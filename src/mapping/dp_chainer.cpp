// SPDX-License-Identifier: MIT

#include "mapping/dp_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::mapping {

namespace {

/* Anchor merging -- merge overlapping PathAnchors on the same diagonal. */

bool anchors_can_merge(const PathAnchor& a, const PathAnchor& b) {
  auto a_diag = static_cast<std::int64_t>(a.ref_coord) - static_cast<std::int64_t>(a.query_pos);
  auto b_diag = static_cast<std::int64_t>(b.ref_coord) - static_cast<std::int64_t>(b.query_pos);
  if (a_diag != b_diag) return false;
  return b.query_pos <= a.query_pos + a.length;
}

void anchors_merge_into(PathAnchor& a, const PathAnchor& b) {
  auto merged_end =
      std::max(static_cast<std::uint32_t>(a.query_pos + a.length),
               static_cast<std::uint32_t>(b.query_pos + b.length));
  a.length = static_cast<std::uint16_t>(std::min(merged_end - a.query_pos, 0xFFFFu));
}

// Merge overlapping anchors within a single path's anchor list.
// Anchors must already be on the same path (caller guarantees via grouping).
void merge_path_anchors(std::vector<PathAnchor>& anchors) {
  if (anchors.size() < 2) return;

  std::sort(anchors.begin(), anchors.end(), [](const PathAnchor& a, const PathAnchor& b) {
    return a.ref_coord != b.ref_coord ? a.ref_coord < b.ref_coord : a.query_pos < b.query_pos;
  });

  std::vector<PathAnchor> merged;
  merged.reserve(anchors.size());

  PathAnchor current = anchors[0];
  for (std::size_t i = 1; i < anchors.size(); ++i) {
    if (anchors_can_merge(current, anchors[i])) {
      anchors_merge_into(current, anchors[i]);
    } else {
      merged.push_back(current);
      current = anchors[i];
    }
  }
  merged.push_back(current);
  anchors = std::move(merged);
}

/* Comparator for sorting anchors by (ref_coord, query_pos). */

struct AnchorComparator {
  bool operator()(const PathAnchor& a, const PathAnchor& b) const {
    if (a.ref_coord != b.ref_coord) return a.ref_coord < b.ref_coord;
    return a.query_pos < b.query_pos;
  }
};

}  // namespace

DPChainer::DPChainer(DPChainerConfig config,
                     const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     const std::vector<std::size_t>& path_lengths)
    : config_(std::move(config)), coords_(coords), path_lengths_(path_lengths) {}

PathAnchorGroups DPChainer::expand(const std::vector<NodeAnchor>& hits) const {
  PathAnchorGroups groups;

  for (std::size_t hit_idx = 0; hit_idx < hits.size(); ++hit_idx) {
    const auto& hit = hits[hit_idx];
    const std::size_t node_id = hit.target.node_id;
    if (node_id >= coords_.size()) continue;

    const auto& node_coords = coords_[node_id];
    if (node_coords.empty()) continue;

    // Clamp span to uint16_t
    auto span = static_cast<std::uint16_t>(std::min(hit.span, static_cast<std::size_t>(0xFFFF)));

    for (const auto& coord : node_coords) {
      std::int64_t ref = coord.ref_coord + static_cast<std::int64_t>(hit.target.offset);
      std::int64_t anchor_end = ref + static_cast<std::int64_t>(span);

      // Skip anchors that extend past path boundary
      if (coord.path_id < path_lengths_.size()) {
        auto path_len = static_cast<std::int64_t>(path_lengths_[coord.path_id]);
        if (anchor_end > path_len || ref < 0) continue;
      }

      // Ensure groups vector is large enough
      if (coord.path_id >= groups.size()) {
        groups.resize(coord.path_id + 1);
      }

      groups[coord.path_id].push_back(PathAnchor{
          .ref_coord = static_cast<std::uint32_t>(ref),
          .query_pos = static_cast<std::uint32_t>(hit.read_pos),
          .length = span,
          ._pad = 0,
          .src_idx = static_cast<std::uint32_t>(hit_idx),
      });
    }
  }
  return groups;
}

ChainResult DPChainer::chain(const std::vector<NodeAnchor>& hits) const {
  auto groups = expand(hits);

  // Count total expanded anchors
  std::size_t total_anchors = 0;
  for (auto& group : groups) {
    if (config_.merge_anchors) {
      merge_path_anchors(group);
    }
    total_anchors += group.size();
  }

  auto result = chain_grouped(groups, hits);
  result.expanded_anchor_count = total_anchors;
  return result;
}

ChainResult DPChainer::chain_grouped(const PathAnchorGroups& groups,
                                     const std::vector<NodeAnchor>& hits) const {
  // Run DP per path, collect all chains
  std::vector<Chain> all_chains;

  for (std::size_t path_id = 0; path_id < groups.size(); ++path_id) {
    const auto& anchors = groups[path_id];
    if (anchors.empty()) continue;

    auto path_chains = chain_one_path(anchors, path_id, hits);
    for (auto& c : path_chains) {
      all_chains.push_back(std::move(c));
    }
  }

  // Sort all chains by score descending
  std::sort(all_chains.begin(), all_chains.end(),
            [](const Chain& a, const Chain& b) { return a.score > b.score; });

  // Keep top max_chains
  if (all_chains.size() > config_.max_chains) {
    all_chains.resize(config_.max_chains);
  }

  // Reassign chain_ids after ranking
  for (std::size_t i = 0; i < all_chains.size(); ++i) {
    for (auto& anchor : all_chains[i].anchors) {
      anchor.chain_id = i;
    }
  }

  // Build result
  ChainResult result;
  result.chains = std::move(all_chains);

  if (!result.chains.empty()) {
    result.score = result.chains[0].score;
    result.anchors = result.chains[0].anchors;
  }

  return result;
}

std::vector<Chain> DPChainer::chain_one_path(const std::vector<PathAnchor>& anchors,
                                             std::size_t path_id,
                                             const std::vector<NodeAnchor>& hits) const {
  /* Sort anchors by (ref_coord, query_pos) */
  std::vector<PathAnchor> sorted = anchors;
  std::sort(sorted.begin(), sorted.end(), AnchorComparator{});

  const std::size_t n = sorted.size();

  /* DP initialization */
  std::vector<double> dp(n, 0.0);
  std::vector<int> pred(n, -1);

  /* DP loop */
  for (std::size_t i = 0; i < n; ++i) {
    const auto& anchor_i = sorted[i];
    double best_score = anchor_score(anchor_i);
    int best_pred = -1;

    std::size_t num_skipped = 0;
    for (std::size_t j = i; j > 0 && num_skipped < config_.max_skip;) {
      --j;
      const auto& anchor_j = sorted[j];

      // Band exceeded
      if (anchor_i.ref_coord - anchor_j.ref_coord > static_cast<std::uint32_t>(config_.max_dist)) {
        break;
      }

      if (!can_chain(anchor_j, anchor_i)) {
        ++num_skipped;
        continue;
      }

      num_skipped = 0;

      // minimap2-style matching bonus
      const double cost = gap_cost(anchor_j, anchor_i);
      auto dq = static_cast<std::int64_t>(anchor_i.query_pos) -
                static_cast<std::int64_t>(anchor_j.query_pos);
      auto dr = static_cast<std::int64_t>(anchor_i.ref_coord) -
                static_cast<std::int64_t>(anchor_j.ref_coord);
      double match_bonus =
          std::min({anchor_score(anchor_i), static_cast<double>(std::max<std::int64_t>(0, dq)),
                    static_cast<double>(std::max<std::int64_t>(0, dr))});
      double score = dp[j] + match_bonus - cost;

      if (score > best_score) {
        best_score = score;
        best_pred = static_cast<int>(j);
      }
    }

    dp[i] = best_score;
    pred[i] = best_pred;
  }

  /* Multi-chain extraction for this path */
  std::vector<Chain> chains;
  std::vector<bool> used(n, false);

  while (chains.size() < config_.max_chains) {
    // Find best unused endpoint
    std::size_t best_idx = 0;
    double best_dp_score = -std::numeric_limits<double>::infinity();
    bool found = false;

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) continue;
      if (dp[i] > best_dp_score) {
        best_dp_score = dp[i];
        best_idx = i;
        found = true;
      }
    }

    if (!found || best_dp_score < static_cast<double>(config_.min_chain_score)) {
      break;
    }

    // Backtrack
    auto chain_indices = backtrack_chain(pred, best_idx);

    // Mark used
    for (std::size_t idx : chain_indices) {
      used[idx] = true;
    }

    // Build Chain with ChainedAnchors (recover node info via src_idx)
    Chain chain;
    chain.score = best_dp_score;
    chain.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
      const auto& anchor = sorted[idx];
      const auto& src = hits[anchor.src_idx];

      ChainedAnchor ca;
      ca.target.node_id = src.target.node_id;
      ca.target.offset = src.target.offset;
      ca.target.length = anchor.length;
      ca.read_pos = anchor.query_pos;
      ca.score = dp[idx];
      ca.chain_id = 0;  // reassigned later across paths
      ca.path_id = path_id;
      ca.ref_coord = anchor.ref_coord;

      chain.anchors.push_back(ca);
    }

    chains.push_back(std::move(chain));
  }

  return chains;
}

bool DPChainer::can_chain(const PathAnchor& j, const PathAnchor& i) const {
  // Query positions must strictly increase
  if (i.query_pos <= j.query_pos) return false;

  auto delta_ref =
      static_cast<std::int64_t>(i.ref_coord) - static_cast<std::int64_t>(j.ref_coord);
  auto delta_query =
      static_cast<std::int64_t>(i.query_pos) - static_cast<std::int64_t>(j.query_pos);

  // Must be forward in both dimensions
  if (delta_ref < 0 || delta_query < 0) return false;

  if (static_cast<std::size_t>(delta_ref) > config_.max_dist ||
      static_cast<std::size_t>(delta_query) > config_.max_dist) {
    return false;
  }

  // Diagonal deviation
  auto diag_dev = std::abs(delta_ref - delta_query);
  if (static_cast<std::size_t>(diag_dev) > config_.max_diag_dev) return false;

  return true;
}

double DPChainer::gap_cost(const PathAnchor& j, const PathAnchor& i) const {
  auto j_ref_end = static_cast<std::int64_t>(j.ref_coord) + static_cast<std::int64_t>(j.length);
  auto j_query_end = static_cast<std::int64_t>(j.query_pos + j.length);

  auto ref_gap = static_cast<std::int64_t>(i.ref_coord) - j_ref_end;
  auto query_gap = static_cast<std::int64_t>(i.query_pos) - j_query_end;

  auto ref_gap_abs = std::max<std::int64_t>(0, ref_gap);
  auto query_gap_abs = std::max<std::int64_t>(0, query_gap);

  double avg_gap = (ref_gap_abs + query_gap_abs) / 2.0;
  double cost = avg_gap * config_.gap_penalty_factor;

  // Diagonal deviation
  auto delta_ref =
      static_cast<std::int64_t>(i.ref_coord) - static_cast<std::int64_t>(j.ref_coord);
  auto delta_query =
      static_cast<std::int64_t>(i.query_pos) - static_cast<std::int64_t>(j.query_pos);
  double diag_dev = std::abs(delta_ref - delta_query);
  cost += diag_dev * config_.diag_penalty_factor;

  // Overlap penalty
  if (ref_gap < 0 || query_gap < 0) {
    double ref_overlap = std::abs(std::min<std::int64_t>(0, ref_gap));
    double query_overlap = std::abs(std::min<std::int64_t>(0, query_gap));
    double avg_overlap = (ref_overlap + query_overlap) / 2.0;
    cost += avg_overlap * config_.overlap_penalty_factor;
  }

  return cost;
}

double DPChainer::anchor_score(const PathAnchor& anchor) const {
  return static_cast<double>(anchor.length) * config_.anchor_weight;
}

std::vector<std::size_t> DPChainer::backtrack_chain(const std::vector<int>& pred,
                                                    std::size_t best_idx) const {
  std::vector<std::size_t> chain;
  int current = static_cast<int>(best_idx);
  while (current != -1) {
    chain.push_back(static_cast<std::size_t>(current));
    current = pred[current];
  }
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
