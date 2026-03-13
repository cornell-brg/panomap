// SPDX-License-Identifier: MIT

#include "mapping/path_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::mapping {

namespace {

// Fast approximate log2 from minimap2/RawHash2 (lchain.c).
// Only valid for x >= 2.
inline float mg_log2(float x) {
  union {
    float f;
    std::uint32_t i;
  } z = {x};
  float log_2 = ((z.i >> 23) & 255) - 128;
  z.i &= ~(255 << 23);
  z.i += 127 << 23;
  log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
  return log_2;
}

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

PathChainer::PathChainer(PathChainerConfig config,
                     const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     const std::vector<std::size_t>& path_lengths)
    : config_(std::move(config)), coords_(coords), path_lengths_(path_lengths) {}

PathAnchorGroups PathChainer::expand(const std::vector<NodeAnchor>& hits) const {
  PathAnchorGroups groups;

  for (std::size_t hit_idx = 0; hit_idx < hits.size(); ++hit_idx) {
    const auto& hit = hits[hit_idx];
    if (hit.node_id >= coords_.size()) continue;

    const auto& node_coords = coords_[hit.node_id];
    if (node_coords.empty()) continue;

    for (const auto& coord : node_coords) {
      std::int64_t ref = coord.ref_coord + static_cast<std::int64_t>(hit.offset);
      std::int64_t anchor_end = ref + static_cast<std::int64_t>(hit.span);

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
          .length = hit.span,
          ._pad = 0,
          .src_idx = static_cast<std::uint32_t>(hit_idx),
      });
    }
  }
  return groups;
}

ChainResult PathChainer::chain(const std::vector<NodeAnchor>& hits) const {
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

ChainResult PathChainer::chain_grouped(const PathAnchorGroups& groups,
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

std::vector<Chain> PathChainer::chain_one_path(const std::vector<PathAnchor>& anchors,
                                             std::size_t path_id,
                                             const std::vector<NodeAnchor>& hits) const {
  /* Sort anchors by (ref_coord, query_pos) */
  std::vector<PathAnchor> sorted = anchors;
  std::sort(sorted.begin(), sorted.end(), AnchorComparator{});

  const std::size_t n = sorted.size();

  /* DP arrays (RawHash2/minimap2 style) */
  std::vector<std::int32_t> f(n, 0);  // f[i] = best score ending at anchor i
  std::vector<std::int32_t> v(n, 0);  // v[i] = peak score up to anchor i
  std::vector<int> pred(n, -1);
  std::vector<int> t(n, 0);  // t[j] = last i that used p[j] as predecessor

  /* DP loop (matches mg_lchain_dp) */
  std::int32_t mmax_f = 0;
  std::size_t st = 0;
  std::int64_t max_ii = -1;

  for (std::size_t i = 0; i < n; ++i) {
    std::int32_t q_span = static_cast<std::int32_t>(sorted[i].length);
    std::int32_t max_f = q_span;  // standalone anchor score
    int max_j = -1;

    // Advance start pointer: skip anchors on different targets or too far away
    while (st < i && sorted[i].ref_coord > sorted[st].ref_coord + config_.max_dist_ref) ++st;

    std::size_t num_skipped = 0;
    for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
      --j;
      auto sc = compute_score(sorted[j], sorted[i]);
      if (sc == std::numeric_limits<std::int32_t>::min()) {
        ++num_skipped;
        continue;
      }

      sc += f[j];
      if (sc > max_f) {
        max_f = sc;
        max_j = static_cast<int>(j);
        if (num_skipped > 0) --num_skipped;
      } else if (t[j] == static_cast<int>(i)) {
        if (++num_skipped > config_.max_skip) break;
      }
      if (pred[j] >= 0) t[pred[j]] = static_cast<int>(i);
    }

    // Check max_ii (best predecessor within band, may have been skipped by the loop)
    if (max_ii < 0 ||
        sorted[i].ref_coord > sorted[max_ii].ref_coord + config_.max_dist_ref) {
      std::int32_t best = std::numeric_limits<std::int32_t>::min();
      max_ii = -1;
      for (std::size_t j = st; j < i; ++j) {
        if (f[j] > best) {
          best = f[j];
          max_ii = static_cast<std::int64_t>(j);
        }
      }
    }

    if (max_ii >= 0 && static_cast<std::size_t>(max_ii) < st) {
      // max_ii fell out of range, skip
    } else if (max_ii >= 0) {
      auto tmp = compute_score(sorted[max_ii], sorted[i]);
      if (tmp != std::numeric_limits<std::int32_t>::min() && max_f < tmp + f[max_ii]) {
        max_f = tmp + f[max_ii];
        max_j = static_cast<int>(max_ii);
      }
    }

    f[i] = max_f;
    pred[i] = max_j;
    v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
    if (max_ii < 0 ||
        (sorted[i].ref_coord <= sorted[max_ii].ref_coord + config_.max_dist_ref &&
         f[max_ii] < f[i])) {
      max_ii = static_cast<std::int64_t>(i);
    }
    if (mmax_f < max_f) mmax_f = max_f;
  }

  /* Multi-chain extraction (backtrack from highest-scoring endpoints) */
  std::vector<Chain> chains;
  std::vector<bool> used(n, false);

  while (chains.size() < config_.max_chains) {
    // Find best unused endpoint
    std::size_t best_idx = 0;
    std::int32_t best_f = std::numeric_limits<std::int32_t>::min();
    bool found = false;

    for (std::size_t i = 0; i < n; ++i) {
      if (used[i]) continue;
      if (f[i] > best_f) {
        best_f = f[i];
        best_idx = i;
        found = true;
      }
    }

    if (!found || best_f < static_cast<std::int32_t>(config_.min_chain_score)) {
      break;
    }

    // Backtrack
    auto chain_indices = backtrack_chain(pred, best_idx);

    // Skip chains with too few anchors
    if (chain_indices.size() < config_.min_chain_anchors) {
      // Mark as used so we don't retry
      for (std::size_t idx : chain_indices) used[idx] = true;
      continue;
    }

    // Mark used
    for (std::size_t idx : chain_indices) {
      used[idx] = true;
    }

    // Build Chain with ChainedAnchors (recover node info via src_idx)
    Chain chain;
    chain.score = static_cast<double>(best_f);
    chain.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
      const auto& anchor = sorted[idx];
      const auto& src = hits[anchor.src_idx];

      ChainedAnchor ca;
      ca.node_id = src.node_id;
      ca.offset = src.offset;
      ca.length = anchor.length;
      ca.read_pos = anchor.query_pos;
      ca.score = static_cast<double>(f[idx]);
      ca.chain_id = 0;  // reassigned later across paths
      ca.path_id = path_id;
      ca.ref_coord = anchor.ref_coord;

      chain.anchors.push_back(ca);
    }

    chains.push_back(std::move(chain));
  }

  return chains;
}

// RawHash2/minimap2-style pairwise scoring (lchain.c:compute_score).
// Returns INT32_MIN if the pair is unchainable.
std::int32_t PathChainer::compute_score(const PathAnchor& j, const PathAnchor& i) const {
  auto dq = static_cast<std::int32_t>(i.query_pos) - static_cast<std::int32_t>(j.query_pos);
  if (dq <= 0 || dq > static_cast<std::int32_t>(config_.max_dist_query))
    return std::numeric_limits<std::int32_t>::min();

  auto dr = static_cast<std::int32_t>(i.ref_coord) - static_cast<std::int32_t>(j.ref_coord);
  if (dr <= 0 || dr > static_cast<std::int32_t>(config_.max_dist_ref))
    return std::numeric_limits<std::int32_t>::min();

  // Diagonal deviation
  auto dd = dr > dq ? dr - dq : dq - dr;
  if (dd > static_cast<std::int32_t>(config_.bw))
    return std::numeric_limits<std::int32_t>::min();

  // Match score: min of q_span and the shorter gap
  auto dg = dr < dq ? dr : dq;
  std::int32_t q_span = static_cast<std::int32_t>(j.length);
  std::int32_t sc = q_span < dg ? q_span : dg;

  // Penalty for diagonal deviation and gap
  if (dd || dg > q_span) {
    float lin_pen = config_.chn_pen_gap * static_cast<float>(dd)
                  + config_.chn_pen_skip * static_cast<float>(dg);
    float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
    sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
  }

  return sc;
}

std::vector<std::size_t> PathChainer::backtrack_chain(const std::vector<int>& pred,
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

// -- PathChainerConfig CLI integration --

std::vector<cli::Option> PathChainerConfig::cli_options() {
  return {
      {'\0', "chain-max-dist", true,
       "DP chain: max query/ref distance for chaining (default: 2500)"},
      {'\0', "chain-bw", true, "DP chain: max diagonal deviation / bandwidth (default: 500)"},
      {'\0', "chain-pen-gap", true, "DP chain: gap penalty factor (default: 0.8)"},
      {'\0', "chain-pen-skip", true, "DP chain: skip penalty factor (default: 0.0)"},
      {'\0', "chain-min-score", true, "DP chain: minimum chain score (default: 15)"},
      {'\0', "chain-min-anchors", true, "DP chain: minimum anchors per chain (default: 2)"},
      {'\0', "chain-max-chains", true, "DP chain: max chains to extract (default: 10)"},
      {'\0', "chain-max-skip", true, "DP chain: stop after N consecutive skips (default: 25)"},
      {'\0', "chain-merge", true, "DP chain: merge overlapping chains (default: false)"},
  };
}

PathChainerConfig PathChainerConfig::from_parsed(const cli::Parsed& parsed) {
  PathChainerConfig cfg;
  if (parsed.values.count("chain-max-dist")) {
    auto val = std::stoull(parsed.values.at("chain-max-dist"));
    cfg.max_dist_ref = val;
    cfg.max_dist_query = val;
  }
  if (parsed.values.count("chain-bw"))
    cfg.bw = std::stoull(parsed.values.at("chain-bw"));
  if (parsed.values.count("chain-pen-gap"))
    cfg.chn_pen_gap = std::stof(parsed.values.at("chain-pen-gap"));
  if (parsed.values.count("chain-pen-skip"))
    cfg.chn_pen_skip = std::stof(parsed.values.at("chain-pen-skip"));
  if (parsed.values.count("chain-min-score"))
    cfg.min_chain_score = std::stoull(parsed.values.at("chain-min-score"));
  if (parsed.values.count("chain-min-anchors"))
    cfg.min_chain_anchors = std::stoull(parsed.values.at("chain-min-anchors"));
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
