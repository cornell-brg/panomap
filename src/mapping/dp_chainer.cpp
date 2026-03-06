// SPDX-License-Identifier: MIT

#include "mapping/dp_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::mapping {

namespace {

/* Anchor merging -- merge overlapping PathAnchors on the same diagonal. */

bool anchors_can_merge(const PathAnchor& a, const PathAnchor& b) {
  // Must be on same diagonal (diagonal = ref_coord - query_pos)
  const std::int64_t a_diag = a.ref_coord - static_cast<std::int64_t>(a.query_pos);
  const std::int64_t b_diag = b.ref_coord - static_cast<std::int64_t>(b.query_pos);
  if (a_diag != b_diag) return false;

  // On same diagonal, check overlap in query
  return b.query_pos <= a.query_pos + a.length;
}

void anchors_merge_into(PathAnchor& a, const PathAnchor& b) {
  const std::size_t merged_end = std::max(a.query_pos + a.length, b.query_pos + b.length);
  a.length = merged_end - a.query_pos;
}

std::vector<PathAnchor> merge_path_anchors(const std::vector<PathAnchor>& anchors) {
  if (anchors.empty()) return {};

  // Group by path_id
  std::map<std::size_t, std::vector<PathAnchor>> by_path;
  for (const auto& a : anchors) by_path[a.path_id].push_back(a);

  std::vector<PathAnchor> merged;
  merged.reserve(anchors.size());

  for (auto& [path_id, path_anchors] : by_path) {
    // Sort by (ref_coord, query_pos)
    std::sort(
        path_anchors.begin(), path_anchors.end(), [](const PathAnchor& a, const PathAnchor& b) {
          return a.ref_coord != b.ref_coord ? a.ref_coord < b.ref_coord : a.query_pos < b.query_pos;
        });

    PathAnchor current = path_anchors[0];
    for (std::size_t i = 1; i < path_anchors.size(); ++i) {
      if (anchors_can_merge(current, path_anchors[i])) {
        anchors_merge_into(current, path_anchors[i]);
      } else {
        merged.push_back(current);
        current = path_anchors[i];
      }
    }
    merged.push_back(current);
  }
  return merged;
}

/* Chain merging -- merge overlapping chains on the same path. */

// Merge overlapping chains on the same path.
// If c1 and c2 overlap (c2 at higher positions), merge them:
// - Combined anchors = union of anchors from both chains
// - Score uses density-based calculation:
//   - c1_only region: c1's score/bp
//   - overlap region: max(c1's score/bp, c2's score/bp)
//   - c2_only region: c2's score/bp
void merge_overlapping_chains(std::vector<Chain>& chains, double /*anchor_weight*/) {
  if (chains.size() < 2) {
    return;
  }

  bool merged = true;
  while (merged) {
    merged = false;

    // Sort chains by (path_id, ref_start)
    std::sort(chains.begin(), chains.end(), [](const Chain& a, const Chain& b) {
      if (a.anchors.empty() || b.anchors.empty()) {
        return a.anchors.size() > b.anchors.size();
      }
      if (a.anchors.front().path_id != b.anchors.front().path_id) {
        return a.anchors.front().path_id < b.anchors.front().path_id;
      }
      return a.anchors.front().ref_coord < b.anchors.front().ref_coord;
    });

    // Try to merge adjacent chains on same path
    for (std::size_t i = 0; i + 1 < chains.size(); ++i) {
      auto& c1 = chains[i];
      auto& c2 = chains[i + 1];

      if (c1.anchors.empty() || c2.anchors.empty()) {
        continue;
      }

      // Check same path
      if (c1.anchors.front().path_id != c2.anchors.front().path_id) {
        continue;
      }

      // Get ref intervals for both chains
      std::int64_t c1_start = c1.anchors.front().ref_coord;
      std::int64_t c1_end =
          c1.anchors.back().ref_coord + static_cast<std::int64_t>(c1.anchors.back().target.length);
      std::int64_t c2_start = c2.anchors.front().ref_coord;
      std::int64_t c2_end =
          c2.anchors.back().ref_coord + static_cast<std::int64_t>(c2.anchors.back().target.length);

      // Check overlap: c2 starts before c1 ends
      if (c2_start < c1_end) {
        // Calculate chain lengths and score densities
        double c1_length = static_cast<double>(c1_end - c1_start);
        double c2_length = static_cast<double>(c2_end - c2_start);
        double c1_density = (c1_length > 0) ? c1.score / c1_length : 0.0;
        double c2_density = (c2_length > 0) ? c2.score / c2_length : 0.0;

        // Calculate region lengths
        // c1_only: from c1_start to c2_start (before overlap)
        // overlap: from c2_start to c1_end
        // c2_only: from c1_end to c2_end (after overlap)
        double c1_only_length = static_cast<double>(c2_start - c1_start);
        double overlap_length = static_cast<double>(c1_end - c2_start);
        double c2_only_length = static_cast<double>(c2_end - c1_end);

        // Clamp negative lengths to 0 (edge cases)
        c1_only_length = std::max(0.0, c1_only_length);
        overlap_length = std::max(0.0, overlap_length);
        c2_only_length = std::max(0.0, c2_only_length);

        // Merged score: use each region's appropriate density
        // Overlap uses the higher density (better chain's score/bp)
        double merged_score = c1_only_length * c1_density +
                              overlap_length * std::max(c1_density, c2_density) +
                              c2_only_length * c2_density;

        // Combine and sort anchors
        std::vector<ChainedAnchor> merged_anchors;
        merged_anchors.reserve(c1.anchors.size() + c2.anchors.size());
        merged_anchors.insert(merged_anchors.end(), c1.anchors.begin(), c1.anchors.end());
        merged_anchors.insert(merged_anchors.end(), c2.anchors.begin(), c2.anchors.end());

        // Sort by ref_coord
        std::sort(merged_anchors.begin(), merged_anchors.end(),
                  [](const ChainedAnchor& a, const ChainedAnchor& b) {
                    return a.ref_coord < b.ref_coord;
                  });

        // Remove duplicates (same ref_coord and read_pos)
        auto last = std::unique(merged_anchors.begin(), merged_anchors.end(),
                                [](const ChainedAnchor& a, const ChainedAnchor& b) {
                                  return a.ref_coord == b.ref_coord && a.read_pos == b.read_pos;
                                });
        merged_anchors.erase(last, merged_anchors.end());

        // Update c1 with merged result
        c1.anchors = std::move(merged_anchors);
        c1.score = merged_score;

        // Remove c2
        chains.erase(chains.begin() + static_cast<std::ptrdiff_t>(i) + 1);

        merged = true;
        break;  // Restart loop
      }
    }
  }

  // Sort by score descending (highest score = rank 0)
  std::sort(chains.begin(), chains.end(),
            [](const Chain& a, const Chain& b) { return a.score > b.score; });

  // Reassign cluster IDs after merging (now reflects rank)
  for (std::size_t i = 0; i < chains.size(); ++i) {
    for (auto& anchor : chains[i].anchors) {
      anchor.chain_id = i;
    }
  }
}

// Check if a coordinate falls within any interval in the list.
// Extends each interval by the specified margin on both ends to catch
// endpoints that are just past the boundary (reducing redundant chains).
bool isInCoveredInterval(std::int64_t coord,
                         const std::vector<std::pair<std::int64_t, std::int64_t>>& intervals,
                         std::int64_t margin = 0) {
  for (const auto& interval : intervals) {
    if (coord >= (interval.first - margin) && coord <= (interval.second + margin)) {
      return true;
    }
  }
  return false;
}

// Comparator for sorting anchors by (path_id, ref_coord, query_pos).
struct AnchorComparator {
  bool operator()(const PathAnchor& a, const PathAnchor& b) const {
    if (a.path_id != b.path_id) {
      return a.path_id < b.path_id;
    }
    if (a.ref_coord != b.ref_coord) {
      return a.ref_coord < b.ref_coord;
    }
    return a.query_pos < b.query_pos;
  }
};

}  // namespace

DPChainer::DPChainer(DPChainerConfig config,
                     const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     const std::vector<std::size_t>& path_lengths)
    : config_(std::move(config)), coords_(coords), path_lengths_(path_lengths) {}

std::vector<PathAnchor> DPChainer::expand(const std::vector<NodeAnchor>& hits) const {
  std::vector<PathAnchor> anchors;
  anchors.reserve(hits.size() * 2);

  for (const auto& hit : hits) {
    const std::size_t node_id = hit.target.node_id;
    if (node_id >= coords_.size()) continue;

    const auto& node_coords = coords_[node_id];
    if (node_coords.empty()) continue;

    for (const auto& coord : node_coords) {
      std::int64_t ref_coord = coord.ref_coord + static_cast<std::int64_t>(hit.target.offset);
      std::int64_t anchor_end = ref_coord + static_cast<std::int64_t>(hit.span);

      // Skip anchors that extend past path boundary
      if (coord.path_id < path_lengths_.size()) {
        std::int64_t path_len = static_cast<std::int64_t>(path_lengths_[coord.path_id]);
        if (anchor_end > path_len || ref_coord < 0) continue;
      }

      anchors.push_back(PathAnchor{
          .query_pos = hit.read_pos,
          .ref_coord = ref_coord,
          .length = hit.span,
          .path_id = coord.path_id,
          .node_id = node_id,
          .node_offset = hit.target.offset,
      });
    }
  }
  return anchors;
}

ChainResult DPChainer::chain(const std::vector<NodeAnchor>& hits) const {
  auto anchors = expand(hits);

  if (config_.merge_anchors) {
    anchors = merge_path_anchors(anchors);
  }

  return chain_path_anchors(anchors);
}

ChainResult DPChainer::chain_path_anchors(const std::vector<PathAnchor>& anchors) const {
  ChainResult summary;

  if (anchors.empty()) {
    return summary;
  }

  // Anchors already expanded by caller (AnchorExpander)
  // Track count for output (already expanded)
  summary.expanded_anchor_count = anchors.size();

  // Step 2: Sort anchors by (path_id, ref_coord, query_pos)
  // Make a mutable copy since we need to sort
  std::vector<PathAnchor> sorted_anchors = anchors;
  std::sort(sorted_anchors.begin(), sorted_anchors.end(), AnchorComparator{});

  const std::size_t n = sorted_anchors.size();

  // Step 3: DP initialization
  std::vector<double> dp(n, 0.0);
  std::vector<int> pred(n, -1);  // -1 means no predecessor (chain starts here)

  // Step 4: DP loop - compute best score for each anchor
  for (std::size_t i = 0; i < n; ++i) {
    const auto& anchor_i = sorted_anchors[i];
    double best_score = anchor_score(anchor_i);  // Start new chain at i
    int best_pred = -1;

    // Check potential predecessors j, starting from closest (i-1) going backwards.
    // This enables early termination via banding: once ref distance exceeds max_dist,
    // all earlier anchors on the same path are even further away.
    //
    // Anchors are sorted by (path_id, ref_coord), so:
    // - Same path anchors are contiguous
    // - When we hit a different path_id, all earlier anchors are also different path
    //
    // max_skip heuristic: stop after consecutive failed chain attempts (minimap2-style).
    // If we've checked max_skip anchors and none can chain, further anchors are unlikely
    // to be better predecessors.
    std::size_t num_skipped = 0;
    for (std::size_t j = i; j > 0 && num_skipped < config_.max_skip;) {
      --j;
      const auto& anchor_j = sorted_anchors[j];

      // Different path - all earlier anchors are also different path, break.
      // (Cross-path chaining is handled separately in DEV033, currently disabled)
      if (anchor_j.path_id != anchor_i.path_id) {
        break;
      }

      // Band exceeded - earlier same-path anchors are even further, break.
      if (anchor_i.ref_coord - anchor_j.ref_coord > static_cast<std::int64_t>(config_.max_dist)) {
        break;
      }

      // Check if j can chain to i (query order, diagonal constraints)
      if (!can_chain(anchor_j, anchor_i)) {
        ++num_skipped;
        continue;
      }

      // Successful chain attempt - reset skip counter
      num_skipped = 0;

      // Compute score if we extend chain ending at j with anchor i
      const double cost = gap_cost(anchor_j, anchor_i);
      // minimap2-style matching bonus: min(Δq, Δr, anchor_length)
      // Penalizes drift -- only get credit for the shorter advance
      const std::int64_t dq = static_cast<std::int64_t>(anchor_i.query_pos) -
                              static_cast<std::int64_t>(anchor_j.query_pos);
      const std::int64_t dr = anchor_i.ref_coord - anchor_j.ref_coord;
      const double match_bonus =
          std::min({anchor_score(anchor_i), static_cast<double>(std::max<std::int64_t>(0, dq)),
                    static_cast<double>(std::max<std::int64_t>(0, dr))});
      const double score = dp[j] + match_bonus - cost;

      if (score > best_score) {
        best_score = score;
        best_pred = static_cast<int>(j);
      }
    }

    dp[i] = best_score;
    pred[i] = best_pred;
  }

  PIRU_PROFILE_START(true, "dp-chain-extraction");

  // Step 5: Multi-chain extraction
  // Extract up to max_chains chains, each with a unique endpoint.
  // Chains can share predecessors (common prefixes).
  // Skip endpoints that fall within already-covered intervals on the same path
  // to avoid redundant same-path chains.
  std::unordered_set<std::size_t> used_in_chain;
  std::map<std::size_t, std::vector<std::pair<std::int64_t, std::int64_t>>> covered_intervals;
  std::size_t chain_id = 0;

  while (chain_id < config_.max_chains) {
    // Find best endpoint among UNUSED anchors, also skipping endpoints
    // that fall within covered intervals on the same path
    std::size_t best_idx = 0;
    double best_dp_score = -std::numeric_limits<double>::infinity();
    bool found = false;

    for (std::size_t i = 0; i < n; ++i) {
      // Skip if anchor already used in a chain
      if (used_in_chain.find(i) != used_in_chain.end()) {
        continue;
      }

      // Skip if endpoint falls within (or near) a covered interval on this path.
      // Use a margin to filter chains that differ only slightly at the boundary.
      constexpr std::int64_t kCoveredMargin = 200;  // bp margin around covered intervals
      const auto& anchor = sorted_anchors[i];
      auto it = covered_intervals.find(anchor.path_id);
      if (it != covered_intervals.end() &&
          isInCoveredInterval(anchor.ref_coord, it->second, kCoveredMargin)) {
        continue;
      }

      if (dp[i] > best_dp_score) {
        best_dp_score = dp[i];
        best_idx = i;
        found = true;
      }
    }

    // Stop if no valid anchors or best score below threshold
    if (!found || best_dp_score < static_cast<double>(config_.min_chain_score)) {
      break;
    }

    // Backtrack to extract chain (can follow pred to used anchors)
    auto chain_indices = backtrack_chain(pred, best_idx);

    // Mark ALL anchors in this chain as used
    for (std::size_t idx : chain_indices) {
      used_in_chain.insert(idx);
    }

    // Compute chain's reference interval for covered tracking
    const auto& first_anchor = sorted_anchors[chain_indices.front()];
    const auto& last_anchor = sorted_anchors[chain_indices.back()];
    std::int64_t ref_start = first_anchor.ref_coord;
    std::int64_t ref_end = last_anchor.ref_coord + static_cast<std::int64_t>(last_anchor.length);

    // Add this chain's interval to covered intervals for its path
    covered_intervals[first_anchor.path_id].emplace_back(ref_start, ref_end);

    // Convert chain anchors to Chain format
    Chain group;
    group.score = best_dp_score;
    group.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
      const auto& anchor = sorted_anchors[idx];
      ChainedAnchor seed_anchor;
      seed_anchor.target.node_id = anchor.node_id;
      seed_anchor.target.offset = anchor.node_offset;
      seed_anchor.target.length = anchor.length;
      seed_anchor.read_pos = anchor.query_pos;
      seed_anchor.score = dp[idx];
      seed_anchor.chain_id = chain_id;

      // Preserve linear coordinates for path-walk debugging
      seed_anchor.path_id = anchor.path_id;
      seed_anchor.ref_coord = anchor.ref_coord;

      group.anchors.push_back(seed_anchor);
    }

    summary.chains.push_back(std::move(group));
    ++chain_id;
  }

  PIRU_PROFILE_STOP(true, "dp-chain-extraction");

  // Post-processing: merge overlapping chains on same path
  if (config_.merge_chains && summary.chains.size() > 1) {
    merge_overlapping_chains(summary.chains, config_.anchor_weight);
  }

  // Set summary score to best chain score (first chain)
  if (!summary.chains.empty()) {
    summary.score = summary.chains[0].score;

    // Also populate flat anchors list with best chain for backward compatibility
    summary.anchors = summary.chains[0].anchors;
  }

  return summary;
}

bool DPChainer::can_chain(const PathAnchor& j, const PathAnchor& i) const {
  // Check order: j must come before i in sorted order (already guaranteed by DP loop)
  // Query positions must strictly increase - equal positions are alternatives, not a chain
  if (i.query_pos <= j.query_pos) {
    return false;
  }

  // Must be same path (cross-path chaining not supported, see DEV033)
  if (i.path_id != j.path_id) {
    return false;
  }

  // Compute deltas
  const std::int64_t delta_ref = i.ref_coord - j.ref_coord;
  const std::int64_t delta_query =
      static_cast<std::int64_t>(i.query_pos) - static_cast<std::int64_t>(j.query_pos);

  // Check distance constraints
  if (delta_ref < 0 || delta_query < 0) {
    return false;  // Must be forward in both dimensions
  }

  if (static_cast<std::size_t>(delta_ref) > config_.max_dist ||
      static_cast<std::size_t>(delta_query) > config_.max_dist) {
    return false;
  }

  // Check diagonal deviation constraint
  const std::int64_t diag_dev = std::abs(delta_ref - delta_query);
  if (static_cast<std::size_t>(diag_dev) > config_.max_diag_dev) {
    return false;
  }

  return true;
}

double DPChainer::gap_cost(const PathAnchor& j, const PathAnchor& i) const {
  // Compute gap between end of j and start of i
  const std::int64_t j_ref_end = j.ref_coord + static_cast<std::int64_t>(j.length);
  const std::int64_t j_query_end = static_cast<std::int64_t>(j.query_pos + j.length);

  const std::int64_t ref_gap = i.ref_coord - j_ref_end;
  const std::int64_t query_gap = static_cast<std::int64_t>(i.query_pos) - j_query_end;

  // Gap is 0 if there's overlap (negative gap)
  const std::int64_t ref_gap_abs = std::max<std::int64_t>(0, ref_gap);
  const std::int64_t query_gap_abs = std::max<std::int64_t>(0, query_gap);

  // Average gap distance
  const double avg_gap = (ref_gap_abs + query_gap_abs) / 2.0;

  // Distance penalty
  double cost = avg_gap * config_.gap_penalty_factor;

  // Diagonal deviation penalty
  const std::int64_t delta_ref = i.ref_coord - j.ref_coord;
  const std::int64_t delta_query =
      static_cast<std::int64_t>(i.query_pos) - static_cast<std::int64_t>(j.query_pos);
  const double diag_dev = std::abs(delta_ref - delta_query);
  cost += diag_dev * config_.diag_penalty_factor;

  // Overlap penalty (if anchors overlap)
  if (ref_gap < 0 || query_gap < 0) {
    const double ref_overlap = std::abs(std::min<std::int64_t>(0, ref_gap));
    const double query_overlap = std::abs(std::min<std::int64_t>(0, query_gap));
    const double avg_overlap = (ref_overlap + query_overlap) / 2.0;
    cost += avg_overlap * config_.overlap_penalty_factor;
  }

  return cost;
}

double DPChainer::anchor_score(const PathAnchor& anchor) const {
  // Score based on anchor length (coverage)
  return static_cast<double>(anchor.length) * config_.anchor_weight;
}

std::vector<std::size_t> DPChainer::backtrack_chain(const std::vector<int>& pred,
                                                    std::size_t best_idx) const {
  std::vector<std::size_t> chain;

  // Backtrack from best_idx following predecessor pointers
  int current = static_cast<int>(best_idx);
  while (current != -1) {
    chain.push_back(static_cast<std::size_t>(current));
    current = pred[current];
  }

  // Reverse to get forward order (chain starts at beginning of read)
  std::reverse(chain.begin(), chain.end());

  return chain;
}

// -- DPChainerConfig CLI integration --

std::vector<cli::Option> DPChainerConfig::cli_options() {
  return {
      {'\0', "chain-max-dist", true,
       "DP chain: max query/ref distance for chaining (default: 500)"},
      {'\0', "chain-max-diag-dev", true, "DP chain: max diagonal deviation (default: 500)"},
      {'\0', "chain-gap-penalty", true, "DP chain: gap penalty factor (default: 0.02)"},
      {'\0', "chain-diag-penalty", true, "DP chain: diagonal penalty factor (default: 0.05)"},
      {'\0', "chain-overlap-penalty", true, "DP chain: overlap penalty factor (default: 0.90)"},
      {'\0', "chain-anchor-weight", true, "DP chain: anchor weight (default: 1.0)"},
      {'\0', "chain-min-score", true, "DP chain: minimum chain score (default: 12)"},
      {'\0', "chain-max-chains", true, "DP chain: max chains to extract (default: 10)"},
      {'\0', "chain-max-skip", true, "DP chain: stop after N consecutive skips (default: 25)"},
      {'\0', "chain-merge", true, "DP chain: merge overlapping chains (default: false)"},
  };
}

DPChainerConfig DPChainerConfig::from_parsed(const cli::Parsed& parsed) {
  DPChainerConfig cfg;
  if (parsed.values.count("chain-max-dist"))
    cfg.max_dist = std::stoull(parsed.values.at("chain-max-dist"));
  if (parsed.values.count("chain-max-diag-dev"))
    cfg.max_diag_dev = std::stoull(parsed.values.at("chain-max-diag-dev"));
  if (parsed.values.count("chain-gap-penalty"))
    cfg.gap_penalty_factor = std::stod(parsed.values.at("chain-gap-penalty"));
  if (parsed.values.count("chain-diag-penalty"))
    cfg.diag_penalty_factor = std::stod(parsed.values.at("chain-diag-penalty"));
  if (parsed.values.count("chain-overlap-penalty"))
    cfg.overlap_penalty_factor = std::stod(parsed.values.at("chain-overlap-penalty"));
  if (parsed.values.count("chain-anchor-weight"))
    cfg.anchor_weight = std::stod(parsed.values.at("chain-anchor-weight"));
  if (parsed.values.count("chain-min-score"))
    cfg.min_chain_score = std::stoull(parsed.values.at("chain-min-score"));
  if (parsed.values.count("chain-max-chains"))
    cfg.max_chains = std::stoull(parsed.values.at("chain-max-chains"));
  if (parsed.values.count("chain-max-skip"))
    cfg.max_skip = std::stoull(parsed.values.at("chain-max-skip"));
  if (parsed.values.count("chain-merge")) {
    const std::string val = parsed.values.at("chain-merge");
    cfg.merge_chains = (val == "true" || val == "1" || val == "yes");
  }
  return cfg;
}

}  // namespace piru::mapping
