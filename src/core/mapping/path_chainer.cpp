/**
 * path_chainer.cpp
 *
 * Path-space DP chainer with SoA layout and component-aware interval dedup.
 * Falls back to node-walk hash dedup when 1D coords are not available.
 *
 * Related:
 *  - path_chainer.hpp
 *  - chainer.hpp
 *
 * SPDX-License-Identifier: MIT
 */

#include "core/mapping/path_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

#include "core/util/logging.hpp"
#include "core/util/timing.hpp"

namespace panomap::mapping {

namespace {

/* SoA DP buffers -- reused across paths within a single chain() call.
 * Kept as a local variable in chain() for thread safety. */
struct DPBuffers {
  /* SoA anchor data (populated from PathAnchors) */
  std::vector<std::uint32_t> ref_coord;
  std::vector<std::uint32_t> query_pos;
  std::vector<std::uint16_t> length;
  std::vector<std::uint32_t> src_idx;
  /* Precomputed per-anchor */
  std::vector<std::int32_t> q_span;
  /* DP state */
  std::vector<std::int32_t> f, v;
  std::vector<std::int32_t> pred, t;
  /* Sort permutation */
  std::vector<std::uint32_t> order;

  void resize(std::size_t n);
};

void DPBuffers::resize(std::size_t n) {
  ref_coord.resize(n);
  query_pos.resize(n);
  length.resize(n);
  src_idx.resize(n);
  q_span.resize(n);
  f.resize(n);
  v.resize(n);
  pred.resize(n);
  t.resize(n);
  order.resize(n);
}

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

// Compact anchor in linear path coordinate space (16 bytes).
// Produced by expanding NodeAnchors using linearization coordinates.
// Internal to PathChainer -- not part of the public API.
struct PathAnchor {
  std::uint32_t ref_coord{0};  // Linear position along reference path
  std::uint32_t query_pos{0};  // Position in query/read
  std::uint16_t length{0};     // Coverage length (from seed span)
  std::uint16_t _pad{0};       // Reserved
  std::uint32_t src_idx{0};    // Index into NodeAnchor input array
};

// Grouped PathAnchors by path_id. path_id is implicit (the vector index).
using PathAnchorGroups = std::vector<std::vector<PathAnchor>>;

/* Anchor merging -- merge overlapping PathAnchors on the same diagonal. */

bool anchors_can_merge(const PathAnchor& a, const PathAnchor& b) {
  auto a_diag = static_cast<std::int64_t>(a.ref_coord) - static_cast<std::int64_t>(a.query_pos);
  auto b_diag = static_cast<std::int64_t>(b.ref_coord) - static_cast<std::int64_t>(b.query_pos);
  if (a_diag != b_diag) return false;
  return b.query_pos <= a.query_pos + a.length;
}

void anchors_merge_into(PathAnchor& a, const PathAnchor& b) {
  auto merged_end = std::max(static_cast<std::uint32_t>(a.query_pos + a.length),
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

// Apply a permutation in-place to parallel SoA arrays.
// order[i] = source index for position i after sorting.
void apply_permutation(DPBuffers& dp, std::size_t n) {
  // Use a temp copy approach (simple, cache-friendly for small n).
  // For very large n a cycle-based in-place permutation would avoid the copy,
  // but pangenome path groups are typically <10k anchors.
  std::vector<std::uint32_t> tmp_u32(n);
  std::vector<std::uint16_t> tmp_u16(n);

  // ref_coord
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.ref_coord[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.ref_coord.begin());

  // query_pos
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.query_pos[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.query_pos.begin());

  // src_idx
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.src_idx[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.src_idx.begin());

  // length
  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.length[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.length.begin());
}

// Hash the ordered node_id sequence of a chain's anchors.
// Chains traversing the same nodes (on different paths) produce the same hash.
// Used as fallback when 1D coords / component IDs are not available.
std::size_t hash_node_walk(const std::vector<ChainedAnchor>& anchors) {
  std::size_t h = 14695981039346656037ULL;  // FNV-1a offset basis
  std::uint32_t prev_node = std::numeric_limits<std::uint32_t>::max();
  for (const auto& a : anchors) {
    if (a.node_id == prev_node) continue;
    h ^= static_cast<std::size_t>(a.node_id);
    h *= 1099511628211ULL;  // FNV prime
    prev_node = a.node_id;
  }
  return h;
}

// Compute 1D interval [start, end] for a chain from node_1d_coords.
struct Interval1D {
  float start;
  float end;
  std::uint32_t component;
};

Interval1D chain_interval(const Chain& chain, const std::vector<float>& node_1d_coords,
                          const std::vector<std::uint32_t>& component_ids) {
  float lo = std::numeric_limits<float>::max();
  float hi = std::numeric_limits<float>::lowest();
  for (const auto& a : chain.anchors) {
    float pos = node_1d_coords[a.node_id] + static_cast<float>(a.offset);
    lo = std::min(lo, pos);
    hi = std::max(hi, pos + static_cast<float>(a.length));
  }
  return {lo, hi, component_ids[chain.anchors[0].node_id]};
}

// Reciprocal overlap fraction between two intervals.
// Returns min(overlap/len_a, overlap/len_b).
float reciprocal_overlap(const Interval1D& a, const Interval1D& b) {
  float overlap = std::min(a.end, b.end) - std::max(a.start, b.start);
  if (overlap <= 0.0f) return 0.0f;
  float len_a = a.end - a.start;
  float len_b = b.end - b.start;
  if (len_a <= 0.0f || len_b <= 0.0f) return 0.0f;
  return std::min(overlap / len_a, overlap / len_b);
}

constexpr float kDedupOverlapThreshold = 0.8f;

}  // namespace

PathChainer::PathChainer(PathChainerConfig config,
                         const std::vector<std::vector<index::LinearCoordinate>>& coords,
                         const std::vector<std::size_t>& path_lengths,
                         const std::vector<float>* node_1d_coords,
                         const std::vector<std::uint32_t>* component_ids)
    : config_(std::move(config)),
      coords_(coords),
      path_lengths_(path_lengths),
      node_1d_coords_(node_1d_coords),
      component_ids_(component_ids) {}

/* Expand NodeAnchors to PathAnchors grouped by path_id. */
PathAnchorGroups expand(const std::vector<NodeAnchor>& hits,
                        const std::vector<std::vector<index::LinearCoordinate>>& coords,
                        const std::vector<std::size_t>& path_lengths) {
  PathAnchorGroups groups;

  for (std::size_t hit_idx = 0; hit_idx < hits.size(); ++hit_idx) {
    const auto& hit = hits[hit_idx];
    if (hit.node_id >= coords.size()) continue;

    const auto& node_coords = coords[hit.node_id];
    if (node_coords.empty()) continue;

    for (const auto& coord : node_coords) {
      std::int64_t ref = coord.ref_coord + static_cast<std::int64_t>(hit.offset);
      std::int64_t anchor_end = ref + static_cast<std::int64_t>(hit.span);

      // Skip anchors that extend past path boundary
      if (coord.path_id < path_lengths.size()) {
        auto path_len = static_cast<std::int64_t>(path_lengths[coord.path_id]);
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
  auto groups = expand(hits, coords_, path_lengths_);

  /* Count total expanded anchors */
  std::size_t total_anchors = 0;
  for (auto& group : groups) {
    if (config_.merge_anchors) {
      merge_path_anchors(group);
    }
    total_anchors += group.size();
  }

  /* Run DP per path, collect chains for dedup. */
  bool use_interval_dedup = node_1d_coords_ && component_ids_;
  std::unordered_map<std::size_t, Chain> walk_map;          // fallback: node-walk hash dedup
  std::vector<Chain> kept_chains;                            // interval dedup: kept chains
  std::vector<Interval1D> kept_intervals;                    // interval dedup: their 1D intervals
  std::vector<bool> input_used(hits.size(), false);          // track which input hits are chained
  DPBuffers dp;                                              // Reused across paths within this call

  for (std::size_t path_id = 0; path_id < groups.size(); ++path_id) {
    const auto& anchors = groups[path_id];
    if (anchors.empty()) continue;

    const std::size_t n = anchors.size();

    /* 1. Populate SoA from PathAnchors */
    dp.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      dp.ref_coord[i] = anchors[i].ref_coord;
      dp.query_pos[i] = anchors[i].query_pos;
      dp.length[i] = anchors[i].length;
      dp.src_idx[i] = anchors[i].src_idx;
      dp.order[i] = static_cast<std::uint32_t>(i);
    }

    /* 2. Sort via index permutation, then apply to all SoA arrays */
    std::sort(dp.order.begin(), dp.order.begin() + n, [&dp](std::uint32_t a, std::uint32_t b) {
      if (dp.ref_coord[a] != dp.ref_coord[b]) return dp.ref_coord[a] < dp.ref_coord[b];
      return dp.query_pos[a] < dp.query_pos[b];
    });
    apply_permutation(dp, n);

    /* 3. Precompute q_span per anchor */
    const bool has_pore_k = config_.pore_k > 0;
    const auto pore_k_i32 = static_cast<std::int32_t>(config_.pore_k);
    for (std::size_t i = 0; i < n; ++i) {
      std::int32_t raw_len = static_cast<std::int32_t>(dp.length[i]);
      dp.q_span[i] = has_pore_k ? pore_k_i32 + raw_len - 1 : raw_len;
    }

    /* 4. Initialize DP state */
    for (std::size_t i = 0; i < n; ++i) {
      dp.f[i] = 0;
      dp.v[i] = 0;
      dp.pred[i] = -1;
      dp.t[i] = 0;
    }

    /* Cache config values for inner loop */
    const auto max_dist_ref = static_cast<std::uint32_t>(config_.max_dist_ref);
    const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
    const auto bw = static_cast<std::int32_t>(config_.bw);
    const float chn_pen_gap = config_.chn_pen_gap;
    const float chn_pen_skip = config_.chn_pen_skip;

    // keeps them in registers
    const auto* rc = dp.ref_coord.data();
    const auto* qp = dp.query_pos.data();
    const auto* qs = dp.q_span.data();
    auto* f = dp.f.data();
    auto* v = dp.v.data();
    auto* pred = dp.pred.data();
    auto* t = dp.t.data();

    /* 5. DP loop (matches mg_lchain_dp) */
    std::size_t st = 0;
    std::int64_t max_ii = -1;

    for (std::size_t i = 0; i < n; ++i) {
      std::int32_t max_f = qs[i];  // standalone anchor score
      std::int32_t max_j = -1;

      // Advance start pointer
      while (st < i && rc[i] > rc[st] + max_dist_ref) ++st;

      std::size_t num_skipped = 0;
      std::size_t num_iter = 0;
      const std::size_t max_iter = config_.max_iterations;
      for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
        --j;
        if (max_iter > 0 && ++num_iter > max_iter) break;

        /* Inline compute_score(j, i) */
        auto dq = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[j]);
        if (dq <= 0 || dq > max_dist_query) {
          ++num_skipped;
          continue;
        }
        auto dr = static_cast<std::int32_t>(rc[i]) - static_cast<std::int32_t>(rc[j]);
        if (dr <= 0 || dr > static_cast<std::int32_t>(max_dist_ref)) {
          ++num_skipped;
          continue;
        }
        auto dd = dr > dq ? dr - dq : dq - dr;
        if (dd > bw) {
          ++num_skipped;
          continue;
        }
        auto dg = dr < dq ? dr : dq;
        std::int32_t sc = qs[j] < dg ? qs[j] : dg;
        if (dd || dg > qs[j]) {
          float lin_pen =
              chn_pen_gap * static_cast<float>(dd) + chn_pen_skip * static_cast<float>(dg);
          float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
          sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
        }
        /* End inline compute_score */

        sc += f[j];
        if (sc > max_f) {
          max_f = sc;
          max_j = static_cast<std::int32_t>(j);
          if (num_skipped > 0) --num_skipped;
        } else if (t[j] == static_cast<std::int32_t>(i)) {
          if (++num_skipped > config_.max_skip) break;
        }
        if (pred[j] >= 0) t[pred[j]] = static_cast<std::int32_t>(i);
      }

      /* Check max_ii (best predecessor within band) */
      if (max_ii < 0 || rc[i] > rc[max_ii] + max_dist_ref) {
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
        /* Inline compute_score(max_ii, i) */
        auto dq_m = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[max_ii]);
        auto dr_m = static_cast<std::int32_t>(rc[i]) - static_cast<std::int32_t>(rc[max_ii]);
        std::int32_t tmp = std::numeric_limits<std::int32_t>::min();
        if (dq_m > 0 && dq_m <= max_dist_query && dr_m > 0 &&
            dr_m <= static_cast<std::int32_t>(max_dist_ref)) {
          auto dd_m = dr_m > dq_m ? dr_m - dq_m : dq_m - dr_m;
          if (dd_m <= bw) {
            auto dg_m = dr_m < dq_m ? dr_m : dq_m;
            tmp = qs[max_ii] < dg_m ? qs[max_ii] : dg_m;
            if (dd_m || dg_m > qs[max_ii]) {
              float lin_pen =
                  chn_pen_gap * static_cast<float>(dd_m) + chn_pen_skip * static_cast<float>(dg_m);
              float log_pen = dd_m >= 1 ? mg_log2(static_cast<float>(dd_m + 1)) : 0.0f;
              tmp -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
            }
          }
        }
        if (tmp != std::numeric_limits<std::int32_t>::min() && max_f < tmp + f[max_ii]) {
          max_f = tmp + f[max_ii];
          max_j = static_cast<std::int32_t>(max_ii);
        }
      }

      f[i] = max_f;
      pred[i] = max_j;
      v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
      if (max_ii < 0 || (rc[i] <= rc[max_ii] + max_dist_ref && f[max_ii] < f[i])) {
        max_ii = static_cast<std::int64_t>(i);
      }
    }

    /* 6. Multi-chain extraction */
    struct ScoreIdx {
      std::int32_t score;
      std::size_t idx;
    };

    std::vector<ScoreIdx> z;
    z.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      if (f[i] >= static_cast<std::int32_t>(config_.min_chain_score) && v[i] <= f[i]) {
        z.push_back({f[i], i});
      }
    }
    std::sort(z.begin(), z.end(),
              [](const ScoreIdx& a, const ScoreIdx& b) { return a.score > b.score; });

    std::vector<bool> used(n, false);
    std::size_t path_chain_count = 0;
    const std::size_t survivor_limit = config_.max_survivor_chains;  // 0 = unlimited

    for (std::size_t k = 0; k < z.size(); ++k) {
      // Stop if both decision chains and survivor chains are exhausted
      bool need_decision = (path_chain_count < config_.max_chains);
      bool need_survivor = (survivor_limit == 0 || path_chain_count < survivor_limit);
      if (!need_decision && !need_survivor) break;

      if (used[z[k].idx]) continue;

      /* Backtrack, marking used as we go.
       * If we hit an already-used anchor, this chain shares a prefix
       * with a previously extracted chain -- discard it entirely. */
      std::vector<std::size_t> chain_indices;
      bool shared_prefix = false;
      {
        std::int32_t cur = static_cast<std::int32_t>(z[k].idx);
        while (cur != -1) {
          auto ci = static_cast<std::size_t>(cur);
          if (used[ci]) {
            shared_prefix = true;
            break;
          }
          used[ci] = true;
          chain_indices.push_back(ci);
          cur = pred[cur];
        }
        std::reverse(chain_indices.begin(), chain_indices.end());
      }

      /* Discard chains that share predecessors or have too few anchors */
      if (shared_prefix || chain_indices.size() < config_.min_chain_anchors) {
        continue;
      }

      /* Mark input_used for survivors (all valid chains, not just top X) */
      for (std::size_t idx : chain_indices) {
        input_used[dp.src_idx[idx]] = true;
      }

      /* Build Chain object only for decision chains (top max_chains) */
      if (need_decision) {
        Chain chain;
        chain.score = static_cast<double>(z[k].score);
        chain.path_id = path_id;
        chain.anchors.reserve(chain_indices.size());

        for (std::size_t idx : chain_indices) {
          auto si = dp.src_idx[idx];
          const auto& src = hits[si];

          ChainedAnchor ca;
          ca.node_id = src.node_id;
          ca.offset = src.offset;
          ca.length = dp.length[idx];
          ca.read_pos = dp.query_pos[idx];
          ca.ref_coord = dp.ref_coord[idx];

          chain.anchors.push_back(ca);
        }

        /* Dedup: interval overlap (preferred) or node-walk hash (fallback). */
        if (use_interval_dedup) {
          auto iv = chain_interval(chain, *node_1d_coords_, *component_ids_);
          bool is_dup = false;
          for (std::size_t ci = 0; ci < kept_chains.size(); ++ci) {
            if (kept_intervals[ci].component != iv.component) continue;
            if (reciprocal_overlap(kept_intervals[ci], iv) >= kDedupOverlapThreshold) {
              // Same region: keep higher score
              if (chain.score > kept_chains[ci].score) {
                kept_chains[ci] = std::move(chain);
                kept_intervals[ci] = iv;
              }
              is_dup = true;
              break;
            }
          }
          if (!is_dup) {
            kept_chains.push_back(std::move(chain));
            kept_intervals.push_back(iv);
          }
        } else {
          std::size_t walk_hash = hash_node_walk(chain.anchors);
          auto it = walk_map.find(walk_hash);
          if (it == walk_map.end()) {
            walk_map.emplace(walk_hash, std::move(chain));
          } else if (chain.score > it->second.score) {
            it->second = std::move(chain);
          }
        }
      }
      ++path_chain_count;
    }
  }

  /* Collect deduplicated chains, sort by score descending */
  std::vector<Chain> deduped;
  if (use_interval_dedup) {
    deduped = std::move(kept_chains);
  } else {
    deduped.reserve(walk_map.size());
    for (auto& [hash, chain] : walk_map) {
      deduped.push_back(std::move(chain));
    }
  }
  std::sort(deduped.begin(), deduped.end(),
            [](const Chain& a, const Chain& b) { return a.score > b.score; });

  /* Trim: keep chains scoring >= secondary_ratio * primary, up to hard cap */
  if (!deduped.empty()) {
    double min_score = deduped[0].score * config_.secondary_ratio;
    std::size_t keep = 1;  // always keep primary
    while (keep < deduped.size() && keep < config_.max_chains &&
           deduped[keep].score >= min_score) {
      ++keep;
    }
    deduped.resize(keep);
  }

  ChainResult result;
  result.chains = std::move(deduped);
  result.expanded_anchor_count = total_anchors;
  result.used_inputs = std::move(input_used);
  return result;
}

/* PathChainerConfig CLI integration */

std::vector<cli::Option> PathChainerConfig::cli_options() {
  return {
      {'\0', "chain-max-dist", true,
       "DP chain: max query/ref distance for chaining (default: 2500)"},
      {'\0', "chain-bw", true, "DP chain: max diagonal deviation / bandwidth (default: 500)"},
      {'\0', "chain-pen-gap", true, "DP chain: gap penalty factor (default: 0.8)"},
      {'\0', "chain-pen-skip", true, "DP chain: skip penalty factor (default: 0.0)"},
      {'\0', "chain-min-score", true, "DP chain: minimum chain score (default: 15)"},
      {'\0', "chain-min-anchors", true, "DP chain: minimum anchors per chain (default: 2)"},
      {'\0', "chain-secondary-ratio", true, "DP chain: stop extracting when score < ratio * primary (default: 0.7)"},
      {'\0', "chain-max-survivor-chains", true,
       "DP chain: max chains for cross-chunk survivor marking (default: 0 = unlimited)"},
      {'\0', "chain-max-skip", true, "DP chain: stop after N consecutive skips (default: 25)"},
      {'\0', "chain-max-iter", true,
       "DP chain: max predecessors to check per anchor (default: 0 = unlimited)"},
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
  if (parsed.values.count("chain-bw")) cfg.bw = std::stoull(parsed.values.at("chain-bw"));
  if (parsed.values.count("chain-pen-gap"))
    cfg.chn_pen_gap = std::stof(parsed.values.at("chain-pen-gap"));
  if (parsed.values.count("chain-pen-skip"))
    cfg.chn_pen_skip = std::stof(parsed.values.at("chain-pen-skip"));
  if (parsed.values.count("chain-min-score"))
    cfg.min_chain_score = std::stoull(parsed.values.at("chain-min-score"));
  if (parsed.values.count("chain-min-anchors"))
    cfg.min_chain_anchors = std::stoull(parsed.values.at("chain-min-anchors"));
  if (parsed.values.count("chain-secondary-ratio"))
    cfg.secondary_ratio = std::stof(parsed.values.at("chain-secondary-ratio"));
  if (parsed.values.count("chain-max-survivor-chains"))
    cfg.max_survivor_chains = std::stoull(parsed.values.at("chain-max-survivor-chains"));
  if (parsed.values.count("chain-max-skip"))
    cfg.max_skip = std::stoull(parsed.values.at("chain-max-skip"));
  if (parsed.values.count("chain-max-iter"))
    cfg.max_iterations = std::stoull(parsed.values.at("chain-max-iter"));
  if (parsed.values.count("chain-merge")) {
    const std::string val = parsed.values.at("chain-merge");
    cfg.merge_chains = (val == "true" || val == "1" || val == "yes");
  }
  return cfg;
}

}  // namespace panomap::mapping
