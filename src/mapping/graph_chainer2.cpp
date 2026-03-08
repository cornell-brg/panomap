/**
 * graph_chainer2.cpp
 *
 * Per-path sorted iteration with shared DP table. For each anchor i
 * (read_pos order), scans backwards in each path's ref_coord-sorted list
 * using a position pointer -- no binary search needed.
 *
 * Related:
 *  - graph_chainer2.hpp
 *  - graph_chainer.cpp  (v1: binary search predecessor lookup)
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

  /* Precompute path memberships for each anchor */
  struct AnchorPathInfo {
    std::size_t path_id;
    std::int64_t ref_coord;
  };
  std::vector<std::vector<AnchorPathInfo>> anchor_paths(n);
  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& a = hits[order[ii]];
    if (a.node_id >= coords_.size()) continue;
    for (const auto& lc : coords_[a.node_id]) {
      anchor_paths[ii].push_back(
          {lc.path_id, lc.ref_coord + static_cast<std::int64_t>(a.offset)});
    }
  }

  /* Build per-path sorted lists with back-pointer to dp_idx.
   * Each path list is sorted by ref_coord. We also store each entry's
   * position (rank) within its path list so we can jump to it in O(1). */
  std::size_t num_paths = transposed_.size();
  struct PathEntry {
    std::int64_t ref_coord;
    std::uint32_t read_pos;
    std::size_t dp_idx;
  };
  std::vector<std::vector<PathEntry>> path_lists(num_paths);

  for (std::size_t ii = 0; ii < n; ++ii) {
    for (const auto& pi : anchor_paths[ii]) {
      path_lists[pi.path_id].push_back(
          {pi.ref_coord, hits[order[ii]].read_pos, ii});
    }
  }
  for (auto& plist : path_lists) {
    std::sort(plist.begin(), plist.end(),
              [](const PathEntry& a, const PathEntry& b) {
                return a.ref_coord < b.ref_coord;
              });
  }

  /* Build position lookup: for each (dp_idx, path_id) -> rank in path_list.
   * Stored as a parallel structure to anchor_paths. */
  std::vector<std::vector<std::size_t>> anchor_path_ranks(n);
  {
    // Count how many entries we've placed per path (to assign ranks)
    std::vector<std::size_t> path_counters(num_paths, 0);

    // We need to match each anchor_paths entry to its position in the
    // sorted path_list. Since path_lists are sorted by ref_coord, we
    // rebuild the rank mapping by iterating each path_list.
    // path_rank_map[path_id][k] = dp_idx of the k-th entry in path_lists[path_id]
    // We invert this to get dp_idx -> rank.

    // Build dp_idx -> list of (path_id, rank) from the sorted path_lists
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> idx_to_ranks(n);
    for (std::size_t p = 0; p < num_paths; ++p) {
      for (std::size_t k = 0; k < path_lists[p].size(); ++k) {
        idx_to_ranks[path_lists[p][k].dp_idx].emplace_back(p, k);
      }
    }

    // Now fill anchor_path_ranks to match anchor_paths ordering
    for (std::size_t ii = 0; ii < n; ++ii) {
      anchor_path_ranks[ii].resize(anchor_paths[ii].size());
      // idx_to_ranks[ii] has (path_id, rank) pairs
      // anchor_paths[ii] has (path_id, ref_coord) pairs -- same order
      // Match by path_id
      for (std::size_t pi_idx = 0; pi_idx < anchor_paths[ii].size(); ++pi_idx) {
        std::size_t target_path = anchor_paths[ii][pi_idx].path_id;
        for (const auto& [pid, rank] : idx_to_ranks[ii]) {
          if (pid == target_path) {
            anchor_path_ranks[ii][pi_idx] = rank;
            break;
          }
        }
      }
    }
  }

  // Debug counters
  std::size_t dbg_candidates = 0;
  std::size_t dbg_valid = 0;
  std::size_t dbg_improved = 0;

  /* DP -- for each anchor i, scan backwards in each path's sorted list */
  std::vector<DPEntry> dp(n);

  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& anchor_i = hits[order[ii]];
    DPEntry best_entry;
    best_entry.score = anchorScore(anchor_i);

    for (std::size_t pi_idx = 0; pi_idx < anchor_paths[ii].size(); ++pi_idx) {
      const auto& pi = anchor_paths[ii][pi_idx];
      const auto& plist = path_lists[pi.path_id];
      std::size_t my_rank = anchor_path_ranks[ii][pi_idx];

      // Scan backwards from my position in this path's sorted list
      std::size_t num_skipped = 0;
      std::int64_t ref_lo = pi.ref_coord - static_cast<std::int64_t>(max_dist_);

      for (std::size_t k = my_rank; k > 0 && num_skipped < max_skip_; --k) {
        const auto& entry = plist[k - 1];
        if (entry.ref_coord < ref_lo) break;

        ++dbg_candidates;

        // Must be earlier in read_pos
        auto delta_query = static_cast<std::int64_t>(anchor_i.read_pos) -
                           static_cast<std::int64_t>(entry.read_pos);
        if (delta_query <= 0) continue;
        if (static_cast<std::size_t>(delta_query) > max_dist_) { ++num_skipped; continue; }

        auto delta_ref = pi.ref_coord - entry.ref_coord;

        auto diag_dev = std::abs(delta_ref - delta_query);
        if (static_cast<std::size_t>(diag_dev) > max_diag_dev_) { ++num_skipped; continue; }

        num_skipped = 0;
        ++dbg_valid;

        // Gap cost
        const auto& anchor_j = hits[order[entry.dp_idx]];
        auto j_ref_end = entry.ref_coord + static_cast<std::int64_t>(anchor_j.span);
        auto j_query_end = static_cast<std::int64_t>(anchor_j.read_pos + anchor_j.span);

        auto ref_gap = pi.ref_coord - j_ref_end;
        auto query_gap = static_cast<std::int64_t>(anchor_i.read_pos) - j_query_end;

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
        if (pi.path_id == dp[entry.dp_idx].path_id) cost -= 0.1;

        // Matching bonus
        double match_bonus =
            std::min({anchorScore(anchor_i),
                      static_cast<double>(std::max<std::int64_t>(0, delta_query)),
                      static_cast<double>(std::max<std::int64_t>(0, delta_ref))});
        double score = dp[entry.dp_idx].score + match_bonus - cost;

        if (score > best_entry.score) {
          best_entry.score = score;
          best_entry.pred = static_cast<int>(entry.dp_idx);
          best_entry.path_id = pi.path_id;
          best_entry.ref_coord = pi.ref_coord;
          ++dbg_improved;
        }
      }
    }

    // If no predecessor found, pick first available path
    if (best_entry.pred == -1 && !anchor_paths[ii].empty()) {
      best_entry.path_id = anchor_paths[ii][0].path_id;
      best_entry.ref_coord = anchor_paths[ii][0].ref_coord;
    }

    dp[ii] = best_entry;
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
