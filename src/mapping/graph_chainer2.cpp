/**
 * graph_chainer2.cpp
 *
 * Path-major DP with shared table. For each path, iterate its ref_coord-sorted
 * anchor list and scan backwards for predecessors (like PathChainer). The DP
 * table is shared across all paths, so scores from earlier paths feed into
 * later ones -- enabling haplotype hopping implicitly.
 *
 * Related:
 *  - graph_chainer2.hpp
 *  - graph_chainer.cpp  (v1: anchor-major with binary search)
 *
 * SPDX-License-Identifier: MIT
 */

#include "mapping/graph_chainer2.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace piru::mapping {

GraphChainer2::GraphChainer2(const std::vector<std::vector<index::LinearCoordinate>>& coords,
                             const std::vector<std::size_t>& path_lengths, std::size_t max_dist,
                             std::size_t max_diag_dev, std::size_t min_chain_score,
                             std::size_t max_chains, std::size_t max_skip)
    : coords_(coords),
      path_lengths_(path_lengths),
      transposed_(index::buildTransposedCoords(coords)),
      max_dist_(max_dist),
      max_diag_dev_(max_diag_dev),
      min_chain_score_(min_chain_score),
      max_chains_(max_chains),
      max_skip_(max_skip) {}

double GraphChainer2::anchorScore(const NodeAnchor& anchor) const {
  return static_cast<double>(anchor.span);
}

ChainResult GraphChainer2::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult result;
  if (hits.empty()) return result;

  /* Sort by read_pos */
  std::vector<std::size_t> order(hits.size());
  for (std::size_t i = 0; i < hits.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return hits[a].read_pos < hits[b].read_pos; });

  const std::size_t n = order.size();
  result.expanded_anchor_count = n;

  /* Per-path sorted lists, self-contained for the inner loop */
  std::size_t num_paths = transposed_.size();
  struct PathEntry {
    std::int64_t ref_coord;
    std::uint32_t read_pos;
    std::uint32_t dp_idx;
    std::uint16_t span;
  };
  std::vector<std::vector<PathEntry>> path_lists(num_paths);

  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& a = hits[order[ii]];
    if (a.node_id >= coords_.size()) continue;
    for (const auto& lc : coords_[a.node_id]) {
      path_lists[lc.path_id].push_back(
          {lc.ref_coord + static_cast<std::int64_t>(a.offset), a.read_pos,
           static_cast<std::uint32_t>(ii), static_cast<std::uint16_t>(a.span)});
    }
  }
  for (auto& plist : path_lists) {
    std::sort(plist.begin(), plist.end(),
              [](const PathEntry& a, const PathEntry& b) {
                return a.ref_coord < b.ref_coord;
              });
  }

  /* Initialize DP with anchor self-scores */
  std::vector<DPEntry> dp(n);
  for (std::size_t ii = 0; ii < n; ++ii) {
    dp[ii].score = anchorScore(hits[order[ii]]);
  }

  // Assign default path for each anchor (first path it appears on)
  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& a = hits[order[ii]];
    if (a.node_id < coords_.size() && !coords_[a.node_id].empty()) {
      const auto& lc = coords_[a.node_id][0];
      dp[ii].path_id = lc.path_id;
      dp[ii].ref_coord = lc.ref_coord + static_cast<std::int64_t>(a.offset);
    }
  }

  // Debug counters
  std::size_t dbg_candidates = 0;
  std::size_t dbg_valid = 0;
  std::size_t dbg_improved = 0;

  /* DP -- path-major: for each path, iterate its sorted list, scan backwards */
  for (std::size_t p = 0; p < num_paths; ++p) {
    const auto& plist = path_lists[p];
    if (plist.empty()) continue;

    for (std::size_t i = 0; i < plist.size(); ++i) {
      const auto& entry_i = plist[i];
      double self_score = static_cast<double>(entry_i.span);

      std::size_t num_skipped = 0;
      std::int64_t ref_lo = entry_i.ref_coord - static_cast<std::int64_t>(max_dist_);

      for (std::size_t k = i; k > 0 && num_skipped < max_skip_; --k) {
        const auto& entry_j = plist[k - 1];
        if (entry_j.ref_coord < ref_lo) break;

        ++dbg_candidates;

        // Must be earlier in read_pos
        auto delta_query = static_cast<std::int64_t>(entry_i.read_pos) -
                           static_cast<std::int64_t>(entry_j.read_pos);
        if (delta_query <= 0) continue;
        if (static_cast<std::size_t>(delta_query) > max_dist_) { ++num_skipped; continue; }

        auto delta_ref = entry_i.ref_coord - entry_j.ref_coord;

        auto diag_dev = std::abs(delta_ref - delta_query);
        if (static_cast<std::size_t>(diag_dev) > max_diag_dev_) { ++num_skipped; continue; }

        num_skipped = 0;
        ++dbg_valid;

        // Gap cost (all fields from PathEntry -- no hits[] access needed)
        auto j_ref_end = entry_j.ref_coord + static_cast<std::int64_t>(entry_j.span);
        auto j_query_end = static_cast<std::int64_t>(entry_j.read_pos + entry_j.span);

        auto ref_gap = entry_i.ref_coord - j_ref_end;
        auto query_gap = static_cast<std::int64_t>(entry_i.read_pos) - j_query_end;

        double ref_gap_abs = std::max<std::int64_t>(0, ref_gap);
        double query_gap_abs = std::max<std::int64_t>(0, query_gap);
        double cost = (ref_gap_abs + query_gap_abs) / 2.0 * 0.02;
        cost += static_cast<double>(diag_dev) * 0.05;

        if (ref_gap < 0 || query_gap < 0) {
          double ref_overlap = std::abs(std::min<std::int64_t>(0, ref_gap));
          double query_overlap = std::abs(std::min<std::int64_t>(0, query_gap));
          cost += (ref_overlap + query_overlap) / 2.0 * 0.90;
        }

        // Same-path preference
        if (p == dp[entry_j.dp_idx].path_id) cost -= 0.1;

        // Matching bonus
        double match_bonus =
            std::min({self_score,
                      static_cast<double>(std::max<std::int64_t>(0, delta_query)),
                      static_cast<double>(std::max<std::int64_t>(0, delta_ref))});
        double score = dp[entry_j.dp_idx].score + match_bonus - cost;

        if (score > dp[entry_i.dp_idx].score) {
          dp[entry_i.dp_idx].score = score;
          dp[entry_i.dp_idx].pred = static_cast<int>(entry_j.dp_idx);
          dp[entry_i.dp_idx].path_id = p;
          dp[entry_i.dp_idx].ref_coord = entry_i.ref_coord;
          ++dbg_improved;
        }
      }
    }
  }

  // Debug summary
  double best_dp = -1e9;
  for (std::size_t i = 0; i < n; ++i) {
    if (dp[i].score > best_dp) best_dp = dp[i].score;
  }
  std::cerr << "[GraphChainer2] anchors=" << n
            << " candidates=" << dbg_candidates
            << " valid=" << dbg_valid
            << " improved=" << dbg_improved
            << " best=" << best_dp
            << "\n";

  /* Multi-chain extraction with truncate-on-used */
  std::vector<bool> used(n, false);

  while (result.chains.size() < max_chains_) {
    // Find best unused endpoint
    std::size_t best_idx = 0;
    double best_dp_score = -std::numeric_limits<double>::infinity();
    bool found = false;

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) continue;
      if (dp[i].score > best_dp_score) {
        best_dp_score = dp[i].score;
        best_idx = i;
        found = true;
      }
    }

    if (!found || best_dp_score < static_cast<double>(min_chain_score_)) break;

    // Backtrack, stop on used anchor
    std::vector<std::size_t> chain_indices;
    {
      int idx = static_cast<int>(best_idx);
      while (idx >= 0) {
        if (used[static_cast<std::size_t>(idx)]) break;
        chain_indices.push_back(static_cast<std::size_t>(idx));
        idx = dp[static_cast<std::size_t>(idx)].pred;
      }
      std::reverse(chain_indices.begin(), chain_indices.end());
    }

    if (chain_indices.empty()) {
      used[best_idx] = true;
      continue;
    }

    // Compute truncated chain score
    double chain_score = best_dp_score;
    std::size_t first_idx = chain_indices.front();
    if (dp[first_idx].pred >= 0) {
      chain_score -= dp[static_cast<std::size_t>(dp[first_idx].pred)].score;
    }

    // Mark used
    for (std::size_t idx : chain_indices) {
      used[idx] = true;
    }

    if (chain_score < static_cast<double>(min_chain_score_)) continue;

    // Build chain
    Chain chain;
    chain.score = chain_score;
    chain.anchors.reserve(chain_indices.size());

    std::size_t chain_path_id = dp[chain_indices.back()].path_id;

    for (std::size_t idx : chain_indices) {
      const auto& src = hits[order[idx]];

      ChainedAnchor ca;
      ca.node_id = src.node_id;
      ca.offset = src.offset;
      ca.length = src.span;
      ca.read_pos = src.read_pos;
      ca.score = dp[idx].score;
      ca.chain_id = result.chains.size();
      ca.path_id = chain_path_id;
      ca.ref_coord = dp[idx].ref_coord;

      chain.anchors.push_back(ca);
    }

    result.chains.push_back(std::move(chain));
  }

  if (!result.chains.empty()) {
    result.score = result.chains[0].score;
    result.anchors = result.chains[0].anchors;
  }

  return result;
}

}  // namespace piru::mapping
