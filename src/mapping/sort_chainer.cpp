/**
 * sort_chainer.cpp
 *
 * 1D sort-based chainer. Converts NodeAnchors to 1D ref coordinates
 * using pre-computed node positions, then runs a single DP pass.
 *
 * Related:
 *  - sort_chainer.hpp
 *  - path_chainer.cpp (reference DP, same SoA approach)
 *
 * SPDX-License-Identifier: MIT
 */

#include "mapping/sort_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Fast approximate log2 (from minimap2/RawHash2)
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

/* SoA buffers for the DP. Separate arrays for each field so the inner
 * loop touches only the arrays it needs (ref_coord, query_pos, q_span). */
struct DPBuffers {
  /* Anchor data */
  std::vector<std::int64_t> ref_coord;
  std::vector<std::uint32_t> query_pos;
  std::vector<std::uint16_t> span;
  std::vector<std::uint32_t> src_idx;  // index back into original hits[]

  /* Precomputed per-anchor */
  std::vector<std::int32_t> q_span;

  /* DP state */
  std::vector<std::int32_t> f, v, pred, t;

  /* Sort permutation */
  std::vector<std::uint32_t> order;

  void resize(std::size_t n) {
    ref_coord.resize(n);
    query_pos.resize(n);
    span.resize(n);
    src_idx.resize(n);
    q_span.resize(n);
    f.resize(n);
    v.resize(n);
    pred.resize(n);
    t.resize(n);
    order.resize(n);
  }
};

// Apply sort permutation in-place to all SoA arrays.
void apply_permutation(DPBuffers& dp, std::size_t n) {
  std::vector<std::int64_t> tmp_i64(n);
  std::vector<std::uint32_t> tmp_u32(n);
  std::vector<std::uint16_t> tmp_u16(n);

  /* ref_coord */
  for (std::size_t i = 0; i < n; ++i) tmp_i64[i] = dp.ref_coord[dp.order[i]];
  std::copy(tmp_i64.begin(), tmp_i64.begin() + n, dp.ref_coord.begin());

  /* query_pos */
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.query_pos[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.query_pos.begin());

  /* src_idx */
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.src_idx[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.src_idx.begin());

  /* span */
  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.span[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.span.begin());
}

}  // namespace

SortChainer::SortChainer(SortChainerConfig config,
                         const std::vector<double>& node_1d_coords)
    : config_(std::move(config)), node_1d_coords_(node_1d_coords) {}

ChainResult SortChainer::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult result;
  if (hits.empty()) return result;

  /* 1. Convert NodeAnchors to SoA arrays */

  DPBuffers dp;
  dp.resize(hits.size());

  std::size_t na = 0;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (h.node_id >= node_1d_coords_.size()) continue;
    double ref = node_1d_coords_[h.node_id] + static_cast<double>(h.offset);
    dp.ref_coord[na] = static_cast<std::int64_t>(ref);
    dp.query_pos[na] = h.read_pos;
    dp.span[na] = h.span;
    dp.src_idx[na] = static_cast<std::uint32_t>(i);
    dp.order[na] = static_cast<std::uint32_t>(na);
    ++na;
  }

  result.expanded_anchor_count = na;
  if (na < config_.min_chain_anchors) return result;

  /* 2. Sort by (ref_coord, query_pos) via index permutation */

  std::sort(dp.order.begin(), dp.order.begin() + na,
            [&dp](std::uint32_t a, std::uint32_t b) {
              if (dp.ref_coord[a] != dp.ref_coord[b])
                return dp.ref_coord[a] < dp.ref_coord[b];
              return dp.query_pos[a] < dp.query_pos[b];
            });
  apply_permutation(dp, na);

  /* 3. Precompute q_span (scoring span per anchor) */

  const bool has_pore_k = config_.pore_k > 0;
  const auto pore_k_i32 = static_cast<std::int32_t>(config_.pore_k);
  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t raw = static_cast<std::int32_t>(dp.span[i]);
    dp.q_span[i] = has_pore_k ? pore_k_i32 + raw - 1 : raw;
  }

  /* 4. Initialize DP state */

  for (std::size_t i = 0; i < na; ++i) {
    dp.f[i] = 0;
    dp.v[i] = 0;
    dp.pred[i] = -1;
    dp.t[i] = 0;
  }

  /* 5. DP forward pass */

  const auto max_dist_ref = static_cast<std::int64_t>(config_.max_dist_ref);
  const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
  const auto bw = static_cast<std::int32_t>(config_.bw);
  const float chn_pen_gap = config_.chn_pen_gap;
  const float chn_pen_skip = config_.chn_pen_skip;

  // Raw pointers for inner loop (keeps them in registers)
  const auto* rc = dp.ref_coord.data();
  const auto* qp = dp.query_pos.data();
  const auto* qs = dp.q_span.data();
  auto* f = dp.f.data();
  auto* v = dp.v.data();
  auto* pred = dp.pred.data();
  auto* t = dp.t.data();

  std::size_t st = 0;
  std::int64_t max_ii = -1;

  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t max_f = qs[i];
    std::int32_t max_j = -1;

    // Advance start pointer past anchors too far back in ref
    while (st < i && rc[i] > rc[st] + max_dist_ref) ++st;

    /* Scan backwards for predecessor candidates */
    std::size_t num_skipped = 0;
    for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
      --j;

      auto dq = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[j]);
      if (dq <= 0 || dq > max_dist_query) { ++num_skipped; continue; }

      auto dr = static_cast<std::int32_t>(rc[i] - rc[j]);
      if (dr <= 0 || dr > static_cast<std::int32_t>(max_dist_ref)) { ++num_skipped; continue; }

      auto dd = dr > dq ? dr - dq : dq - dr;
      if (dd > bw) { ++num_skipped; continue; }

      /* Compute score for j -> i transition */
      auto dg = dr < dq ? dr : dq;
      std::int32_t sc = qs[j] < dg ? qs[j] : dg;
      if (dd || dg > qs[j]) {
        float lin_pen = chn_pen_gap * static_cast<float>(dd)
                      + chn_pen_skip * static_cast<float>(dg);
        float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
        sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
      }

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
        if (f[j] > best) { best = f[j]; max_ii = static_cast<std::int64_t>(j); }
      }
    }

    if (max_ii >= 0 && static_cast<std::size_t>(max_ii) >= st) {
      auto dq_m = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[max_ii]);
      auto dr_m = static_cast<std::int32_t>(rc[i] - rc[max_ii]);
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

  /* 6. Multi-chain extraction
   * Find-max-and-extract loop: O(max_chains * n) which is better than
   * sorting all candidates when n >> max_chains. */

  std::vector<bool> used(na, false);
  std::vector<bool> input_used(hits.size(), false);

  while (result.chains.size() < config_.max_chains) {
    /* Find best unused endpoint (v[i] <= f[i] filters past-peak anchors) */
    std::size_t best_idx = 0;
    std::int32_t best_score = std::numeric_limits<std::int32_t>::min();
    for (std::size_t i = 0; i < na; ++i) {
      if (!used[i] && f[i] >= static_cast<std::int32_t>(config_.min_chain_score)
          && v[i] <= f[i] && f[i] > best_score) {
        best_score = f[i];
        best_idx = i;
      }
    }
    if (best_score < static_cast<std::int32_t>(config_.min_chain_score)) break;

    /* Backtrack, stop on used anchor (shared prefix discard) */
    std::vector<std::size_t> chain_indices;
    bool shared_prefix = false;
    {
      std::int32_t cur = static_cast<std::int32_t>(best_idx);
      while (cur != -1) {
        auto ci = static_cast<std::size_t>(cur);
        if (used[ci]) { shared_prefix = true; break; }
        used[ci] = true;
        chain_indices.push_back(ci);
        cur = pred[cur];
      }
      std::reverse(chain_indices.begin(), chain_indices.end());
    }

    if (shared_prefix || chain_indices.size() < config_.min_chain_anchors) continue;

    /* Build Chain from SoA indices */
    Chain chain;
    chain.score = static_cast<double>(best_score);
    chain.path_id = 0;
    chain.coord_space = CoordSpace::kCanonical;
    chain.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
      auto si = dp.src_idx[idx];
      const auto& src = hits[si];
      input_used[si] = true;

      ChainedAnchor ca;
      ca.node_id = src.node_id;
      ca.offset = src.offset;
      ca.length = dp.span[idx];
      ca.read_pos = dp.query_pos[idx];
      ca.ref_coord = dp.ref_coord[idx];
      chain.anchors.push_back(ca);
    }

    result.chains.push_back(std::move(chain));
  }

  result.used_inputs = std::move(input_used);
  return result;
}

}  // namespace piru::mapping
