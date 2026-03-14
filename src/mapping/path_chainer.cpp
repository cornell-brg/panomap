// SPDX-License-Identifier: MIT

#include "mapping/path_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::mapping {

// -- DPBuffers --

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

}  // namespace

PathChainer::PathChainer(PathChainerConfig config,
                     const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     const std::vector<std::size_t>& path_lengths)
    : config_(std::move(config)), coords_(coords), path_lengths_(path_lengths) {}

// Expand NodeAnchors to PathAnchors grouped by path_id.
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

  // Count total expanded anchors
  std::size_t total_anchors = 0;
  for (auto& group : groups) {
    if (config_.merge_anchors) {
      merge_path_anchors(group);
    }
    total_anchors += group.size();
  }

  // Run DP per path, collect all chains
  std::vector<Chain> all_chains;

  for (std::size_t path_id = 0; path_id < groups.size(); ++path_id) {
    const auto& anchors = groups[path_id];
    if (anchors.empty()) continue;

    const std::size_t n = anchors.size();

    /* 1. Populate SoA from PathAnchors */
    dp_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      dp_.ref_coord[i] = anchors[i].ref_coord;
      dp_.query_pos[i] = anchors[i].query_pos;
      dp_.length[i] = anchors[i].length;
      dp_.src_idx[i] = anchors[i].src_idx;
      dp_.order[i] = static_cast<std::uint32_t>(i);
    }

    /* 2. Sort via index permutation, then apply to all SoA arrays */
    std::sort(dp_.order.begin(), dp_.order.begin() + n,
              [this](std::uint32_t a, std::uint32_t b) {
                if (dp_.ref_coord[a] != dp_.ref_coord[b])
                  return dp_.ref_coord[a] < dp_.ref_coord[b];
                return dp_.query_pos[a] < dp_.query_pos[b];
              });
    apply_permutation(dp_, n);

    /* 3. Precompute q_span per anchor */
    const bool has_pore_k = config_.pore_k > 0;
    const auto pore_k_i32 = static_cast<std::int32_t>(config_.pore_k);
    for (std::size_t i = 0; i < n; ++i) {
      std::int32_t raw_len = static_cast<std::int32_t>(dp_.length[i]);
      dp_.q_span[i] = has_pore_k ? pore_k_i32 + raw_len - 1 : raw_len;
    }

    /* 4. Initialize DP state */
    for (std::size_t i = 0; i < n; ++i) {
      dp_.f[i] = 0;
      dp_.v[i] = 0;
      dp_.pred[i] = -1;
      dp_.t[i] = 0;
    }

    // Cache config values for inner loop
    const auto max_dist_ref = static_cast<std::uint32_t>(config_.max_dist_ref);
    const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
    const auto bw = static_cast<std::int32_t>(config_.bw);
    const float chn_pen_gap = config_.chn_pen_gap;
    const float chn_pen_skip = config_.chn_pen_skip;

    // Local pointers for SoA arrays (keeps them in registers)
    const auto* rc = dp_.ref_coord.data();
    const auto* qp = dp_.query_pos.data();
    const auto* qs = dp_.q_span.data();
    auto* f = dp_.f.data();
    auto* v = dp_.v.data();
    auto* pred = dp_.pred.data();
    auto* t = dp_.t.data();

    /* 5. DP loop (matches mg_lchain_dp) */
    std::size_t st = 0;
    std::int64_t max_ii = -1;

    for (std::size_t i = 0; i < n; ++i) {
      std::int32_t max_f = qs[i];  // standalone anchor score
      std::int32_t max_j = -1;

      // Advance start pointer
      while (st < i && rc[i] > rc[st] + max_dist_ref) ++st;

      std::size_t num_skipped = 0;
      for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
        --j;

        // Inline compute_score(j, i)
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
          float lin_pen = chn_pen_gap * static_cast<float>(dd)
                        + chn_pen_skip * static_cast<float>(dg);
          float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
          sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
        }
        // End inline compute_score

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

      // Check max_ii (best predecessor within band)
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
        // Inline compute_score(max_ii, i)
        auto dq_m = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[max_ii]);
        auto dr_m = static_cast<std::int32_t>(rc[i]) - static_cast<std::int32_t>(rc[max_ii]);
        std::int32_t tmp = std::numeric_limits<std::int32_t>::min();
        if (dq_m > 0 && dq_m <= max_dist_query &&
            dr_m > 0 && dr_m <= static_cast<std::int32_t>(max_dist_ref)) {
          auto dd_m = dr_m > dq_m ? dr_m - dq_m : dq_m - dr_m;
          if (dd_m <= bw) {
            auto dg_m = dr_m < dq_m ? dr_m : dq_m;
            tmp = qs[max_ii] < dg_m ? qs[max_ii] : dg_m;
            if (dd_m || dg_m > qs[max_ii]) {
              float lin_pen = chn_pen_gap * static_cast<float>(dd_m)
                            + chn_pen_skip * static_cast<float>(dg_m);
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
      if (max_ii < 0 ||
          (rc[i] <= rc[max_ii] + max_dist_ref && f[max_ii] < f[i])) {
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

    for (std::size_t k = 0; k < z.size() && path_chain_count < config_.max_chains; ++k) {
      if (used[z[k].idx]) continue;

      // Backtrack, marking used as we go.
      // If we hit an already-used anchor, this chain shares a prefix
      // with a previously extracted chain -- discard it entirely.
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

      // Discard chains that share predecessors with an existing chain,
      // or have too few anchors
      if (shared_prefix || chain_indices.size() < config_.min_chain_anchors) {
        continue;
      }

      // Build Chain: recover graph-space fields via src_idx
      Chain chain;
      chain.score = static_cast<double>(z[k].score);
      chain.path_id = path_id;
      chain.anchors.reserve(chain_indices.size());

      for (std::size_t idx : chain_indices) {
        const auto& src = hits[dp_.src_idx[idx]];

        ChainedAnchor ca;
        ca.node_id = src.node_id;
        ca.offset = src.offset;
        ca.length = dp_.length[idx];
        ca.read_pos = dp_.query_pos[idx];
        ca.ref_coord = dp_.ref_coord[idx];

        chain.anchors.push_back(ca);
      }

      all_chains.push_back(std::move(chain));
      ++path_chain_count;
    }
  }

  // Sort all chains by score descending
  std::sort(all_chains.begin(), all_chains.end(),
            [](const Chain& a, const Chain& b) { return a.score > b.score; });

  // Keep top max_chains
  if (all_chains.size() > config_.max_chains) {
    all_chains.resize(config_.max_chains);
  }

  ChainResult result;
  result.chains = std::move(all_chains);
  result.expanded_anchor_count = total_anchors;
  return result;
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
