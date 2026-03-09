// SPDX-License-Identifier: MIT

#include "mapping/graph_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace piru::mapping {

GraphChainer::GraphChainer(const std::vector<std::vector<index::LinearCoordinate>>& coords,
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

GraphChainer::Transition GraphChainer::bestTransition(const NodeAnchor& j, const NodeAnchor& i,
                                                      std::size_t prev_path_id) const {
  Transition best;
  best.valid = false;
  double best_cost = std::numeric_limits<double>::max();

  if (j.node_id >= coords_.size()) return best;

  const auto& j_coords = coords_[j.node_id];

  for (const auto& lc_j : j_coords) {
    std::int64_t ref_j = lc_j.ref_coord + static_cast<std::int64_t>(j.offset);

    // Look up i's ref_coords on the same path
    if (lc_j.path_id >= transposed_.size()) continue;
    if (i.node_id >= transposed_[lc_j.path_id].size()) continue;

    const auto& i_refs = transposed_[lc_j.path_id][i.node_id];
    if (i_refs.empty()) continue;

    for (std::int64_t ref_coord_i : i_refs) {
      std::int64_t ref_i = ref_coord_i + static_cast<std::int64_t>(i.offset);

      // Ref distance band
      auto delta_ref = ref_i - ref_j;
      if (delta_ref < 0) continue;  // must be forward
      if (static_cast<std::size_t>(delta_ref) > max_dist_) continue;

      // Query distance (already checked by caller, but compute for gap cost)
      auto delta_query = static_cast<std::int64_t>(i.read_pos) - static_cast<std::int64_t>(j.read_pos);
      if (delta_query <= 0) continue;

      // Diagonal deviation
      auto diag_dev = std::abs(delta_ref - delta_query);
      if (static_cast<std::size_t>(diag_dev) > max_diag_dev_) continue;

      // Compute gap cost (same formula as PathChainer)
      auto j_ref_end = ref_j + static_cast<std::int64_t>(j.span);
      auto j_query_end = static_cast<std::int64_t>(j.read_pos + j.span);

      auto ref_gap = ref_i - j_ref_end;
      auto query_gap = static_cast<std::int64_t>(i.read_pos) - j_query_end;

      double ref_gap_abs = std::max<std::int64_t>(0, ref_gap);
      double query_gap_abs = std::max<std::int64_t>(0, query_gap);
      double cost = (ref_gap_abs + query_gap_abs) / 2.0 * 0.02;  // gap_penalty_factor

      cost += static_cast<double>(diag_dev) * 0.05;  // diag_penalty_factor

      // Overlap penalty
      if (ref_gap < 0 || query_gap < 0) {
        double ref_overlap = std::abs(std::min<std::int64_t>(0, ref_gap));
        double query_overlap = std::abs(std::min<std::int64_t>(0, query_gap));
        cost += (ref_overlap + query_overlap) / 2.0 * 0.90;  // overlap_penalty_factor
      }

      // Same-path preference: slight bonus for staying on same path
      if (lc_j.path_id == prev_path_id) {
        cost -= 0.1;
      }

      if (cost < best_cost) {
        best_cost = cost;
        best.cost = cost;
        best.path_id = lc_j.path_id;
        best.ref_j = ref_j;
        best.ref_i = ref_i;
        best.valid = true;
      }
    }
  }

  return best;
}

double GraphChainer::anchorScore(const NodeAnchor& anchor) const {
  return static_cast<double>(anchor.span);
}

std::vector<std::size_t> GraphChainer::backtrack(const std::vector<DPEntry>& dp,
                                                  std::size_t best_idx) const {
  std::vector<std::size_t> indices;
  int idx = static_cast<int>(best_idx);
  while (idx >= 0) {
    indices.push_back(static_cast<std::size_t>(idx));
    idx = dp[static_cast<std::size_t>(idx)].pred;
  }
  std::reverse(indices.begin(), indices.end());
  return indices;
}

ChainResult GraphChainer::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult result;
  if (hits.empty()) return result;

  /* Sort by read_pos */
  std::vector<std::size_t> order(hits.size());
  for (std::size_t i = 0; i < hits.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return hits[a].read_pos < hits[b].read_pos; });

  const std::size_t n = order.size();
  result.expanded_anchor_count = n;

  // Per-path anchor lists sorted by ref_coord, self-contained for the inner
  // loop (no random access back into hits[] needed).
  std::size_t num_paths = transposed_.size();
  struct PathEntry {
    std::int64_t ref_coord;
    std::uint32_t read_pos;
    std::uint32_t dp_idx;
    std::uint16_t span;
    bool operator<(const PathEntry& o) const { return ref_coord < o.ref_coord; }
  };
  std::vector<std::vector<PathEntry>> path_lists(num_paths);

  // Also record which paths each anchor belongs to (for the outer DP loop)
  struct AnchorPathInfo {
    std::uint32_t path_id;
    std::int64_t ref_coord;
  };
  std::vector<std::vector<AnchorPathInfo>> anchor_paths(n);

  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& a = hits[order[ii]];
    if (a.node_id >= coords_.size()) continue;
    for (const auto& lc : coords_[a.node_id]) {
      auto ref = lc.ref_coord + static_cast<std::int64_t>(a.offset);
      path_lists[lc.path_id].push_back(
          {ref, a.read_pos, static_cast<std::uint32_t>(ii),
           static_cast<std::uint16_t>(a.span)});
      anchor_paths[ii].push_back(
          {static_cast<std::uint32_t>(lc.path_id), ref});
    }
  }
  for (auto& plist : path_lists) {
    std::sort(plist.begin(), plist.end());
  }

  // Debug counters
  std::size_t dbg_candidates = 0;
  std::size_t dbg_valid = 0;
  std::size_t dbg_improved = 0;

  /* DP -- path-aware predecessor lookup with binary search */
  std::vector<DPEntry> dp(n);

  for (std::size_t ii = 0; ii < n; ++ii) {
    const auto& anchor_i = hits[order[ii]];
    DPEntry best_entry;
    best_entry.score = anchorScore(anchor_i);

    // For each path that anchor i is on, binary search for predecessors
    for (const auto& pi : anchor_paths[ii]) {
      const auto& plist = path_lists[pi.path_id];
      if (plist.empty()) continue;

      std::int64_t ref_lo = pi.ref_coord - static_cast<std::int64_t>(max_dist_);

      // Binary search for first entry >= ref_lo
      PathEntry search_key{ref_lo, 0, 0, 0};
      auto it_lo = std::lower_bound(plist.begin(), plist.end(), search_key);

      // Scan from lower bound to ref_i
      std::size_t num_skipped = 0;
      for (auto it = it_lo; it != plist.end() && num_skipped < max_skip_; ++it) {
        if (it->ref_coord >= pi.ref_coord) break;

        ++dbg_candidates;

        // Must be a predecessor (earlier in read_pos order)
        auto delta_query = static_cast<std::int64_t>(anchor_i.read_pos) -
                           static_cast<std::int64_t>(it->read_pos);
        if (delta_query <= 0) continue;
        if (static_cast<std::size_t>(delta_query) > max_dist_) { ++num_skipped; continue; }

        auto delta_ref = pi.ref_coord - it->ref_coord;

        auto diag_dev = std::abs(delta_ref - delta_query);
        if (static_cast<std::size_t>(diag_dev) > max_diag_dev_) { ++num_skipped; continue; }

        num_skipped = 0;
        ++dbg_valid;

        // Gap cost (all fields from PathEntry -- no hits[] access needed)
        auto j_ref_end = it->ref_coord + static_cast<std::int64_t>(it->span);
        auto j_query_end = static_cast<std::int64_t>(it->read_pos + it->span);

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
        if (pi.path_id == dp[it->dp_idx].path_id) cost -= 0.1;

        // Matching bonus
        double match_bonus =
            std::min({anchorScore(anchor_i),
                      static_cast<double>(std::max<std::int64_t>(0, delta_query)),
                      static_cast<double>(std::max<std::int64_t>(0, delta_ref))});
        double score = dp[it->dp_idx].score + match_bonus - cost;

        if (score > best_entry.score) {
          best_entry.score = score;
          best_entry.pred = static_cast<int>(it->dp_idx);
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
  std::cerr << "[GraphChainer] anchors=" << n
            << " candidates=" << dbg_candidates
            << " valid=" << dbg_valid
            << " improved=" << dbg_improved
            << " best=" << best_dp
            << "\n";

  /* Multi-chain extraction */
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

    if (!found || best_dp_score < static_cast<double>(min_chain_score_)) {
      break;
    }

    // Backtrack, but stop if we hit a used anchor
    std::vector<std::size_t> chain_indices;
    {
      int idx = static_cast<int>(best_idx);
      while (idx >= 0) {
        if (used[static_cast<std::size_t>(idx)]) break;  // hit previously used anchor
        chain_indices.push_back(static_cast<std::size_t>(idx));
        idx = dp[static_cast<std::size_t>(idx)].pred;
      }
      std::reverse(chain_indices.begin(), chain_indices.end());
    }

    if (chain_indices.empty()) {
      used[best_idx] = true;
      continue;
    }

    // Compute truncated chain score: endpoint score minus the score at
    // the truncation point's predecessor (if truncated)
    double chain_score = best_dp_score;
    std::size_t first_idx = chain_indices.front();
    if (dp[first_idx].pred >= 0) {
      // Truncated -- subtract the score we cut off
      chain_score -= dp[static_cast<std::size_t>(dp[first_idx].pred)].score;
    }

    // Mark chain anchors as used
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

  // Set top-level result from best chain
  if (!result.chains.empty()) {
    result.score = result.chains[0].score;
    result.anchors = result.chains[0].anchors;
  }

  return result;
}

}  // namespace piru::mapping
