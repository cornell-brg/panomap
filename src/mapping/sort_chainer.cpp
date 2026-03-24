/**
 * sort_chainer.cpp
 *
 * 1D sort-based chainer. Converts NodeAnchors to 1D ref coordinates
 * using pre-computed node positions, then runs a single DP pass.
 * DP logic reused from PathChainer (mg_lchain_dp style).
 *
 * Related:
 *  - sort_chainer.hpp
 *  - path_chainer.cpp (reference DP)
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

// Fast approximate log2 (from PathChainer / minimap2)
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

}  // namespace

SortChainer::SortChainer(SortChainerConfig config,
                         const std::vector<double>& node_1d_coords)
    : config_(std::move(config)), node_1d_coords_(node_1d_coords) {}

ChainResult SortChainer::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult result;
  if (hits.empty()) return result;

  const std::size_t n = hits.size();
  result.expanded_anchor_count = n;

  /* 1. Convert NodeAnchors to (1d_ref_coord, query_pos) and sort by ref */
  struct SortAnchor {
    std::int64_t ref_coord;
    std::uint32_t query_pos;
    std::uint16_t span;
    std::uint32_t src_idx;  // index into original hits[]
  };

  std::vector<SortAnchor> anchors;
  anchors.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& h = hits[i];
    if (h.node_id >= node_1d_coords_.size()) continue;
    double ref = node_1d_coords_[h.node_id] + static_cast<double>(h.offset);
    anchors.push_back({
        static_cast<std::int64_t>(ref),
        h.read_pos,
        h.span,
        static_cast<std::uint32_t>(i),
    });
  }

  if (anchors.size() < config_.min_chain_anchors) return result;

  // Sort by ref_coord, then query_pos
  std::sort(anchors.begin(), anchors.end(), [](const SortAnchor& a, const SortAnchor& b) {
    return a.ref_coord != b.ref_coord ? a.ref_coord < b.ref_coord : a.query_pos < b.query_pos;
  });

  const std::size_t na = anchors.size();

  /* 2. DP (same as PathChainer, single pass) */
  std::vector<std::int32_t> f(na), v(na), pred(na), t(na);

  const auto max_dist_ref = static_cast<std::uint32_t>(config_.max_dist_ref);
  const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
  const auto bw = static_cast<std::int32_t>(config_.bw);
  const float chn_pen_gap = config_.chn_pen_gap;
  const float chn_pen_skip = config_.chn_pen_skip;
  const bool has_pore_k = config_.pore_k > 0;
  const auto pore_k_i32 = static_cast<std::int32_t>(config_.pore_k);

  std::size_t st = 0;
  std::int64_t max_ii = -1;

  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t q_span = has_pore_k
        ? pore_k_i32 + static_cast<std::int32_t>(anchors[i].span) - 1
        : static_cast<std::int32_t>(anchors[i].span);

    std::int32_t max_f = q_span;
    std::int32_t max_j = -1;

    while (st < i && anchors[i].ref_coord > anchors[st].ref_coord + max_dist_ref) ++st;

    std::size_t num_skipped = 0;
    for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
      --j;

      auto dq = static_cast<std::int32_t>(anchors[i].query_pos) -
                static_cast<std::int32_t>(anchors[j].query_pos);
      if (dq <= 0 || dq > max_dist_query) { ++num_skipped; continue; }

      auto dr = static_cast<std::int32_t>(anchors[i].ref_coord - anchors[j].ref_coord);
      if (dr <= 0 || dr > static_cast<std::int32_t>(max_dist_ref)) { ++num_skipped; continue; }

      auto dd = dr > dq ? dr - dq : dq - dr;
      if (dd > bw) { ++num_skipped; continue; }

      auto dg = dr < dq ? dr : dq;
      std::int32_t q_span_j = has_pore_k
          ? pore_k_i32 + static_cast<std::int32_t>(anchors[j].span) - 1
          : static_cast<std::int32_t>(anchors[j].span);
      std::int32_t sc = q_span_j < dg ? q_span_j : dg;
      if (dd || dg > q_span_j) {
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

    // Check max_ii
    if (max_ii < 0 || anchors[i].ref_coord > anchors[max_ii].ref_coord + max_dist_ref) {
      std::int32_t best = std::numeric_limits<std::int32_t>::min();
      max_ii = -1;
      for (std::size_t j = st; j < i; ++j) {
        if (f[j] > best) { best = f[j]; max_ii = static_cast<std::int64_t>(j); }
      }
    }

    if (max_ii >= 0 && static_cast<std::size_t>(max_ii) >= st) {
      auto dq_m = static_cast<std::int32_t>(anchors[i].query_pos) -
                  static_cast<std::int32_t>(anchors[max_ii].query_pos);
      auto dr_m = static_cast<std::int32_t>(anchors[i].ref_coord - anchors[max_ii].ref_coord);
      std::int32_t tmp = std::numeric_limits<std::int32_t>::min();
      if (dq_m > 0 && dq_m <= max_dist_query &&
          dr_m > 0 && dr_m <= static_cast<std::int32_t>(max_dist_ref)) {
        auto dd_m = dr_m > dq_m ? dr_m - dq_m : dq_m - dr_m;
        if (dd_m <= bw) {
          std::int32_t q_span_m = has_pore_k
              ? pore_k_i32 + static_cast<std::int32_t>(anchors[max_ii].span) - 1
              : static_cast<std::int32_t>(anchors[max_ii].span);
          auto dg_m = dr_m < dq_m ? dr_m : dq_m;
          tmp = q_span_m < dg_m ? q_span_m : dg_m;
          if (dd_m || dg_m > q_span_m) {
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
        (anchors[i].ref_coord <= anchors[max_ii].ref_coord + max_dist_ref && f[max_ii] < f[i])) {
      max_ii = static_cast<std::int64_t>(i);
    }
  }

  /* 3. Multi-chain extraction */
  struct ScoreIdx {
    std::int32_t score;
    std::size_t idx;
  };

  std::vector<ScoreIdx> z;
  z.reserve(na);
  for (std::size_t i = 0; i < na; ++i) {
    if (f[i] >= static_cast<std::int32_t>(config_.min_chain_score) && v[i] <= f[i]) {
      z.push_back({f[i], i});
    }
  }
  std::sort(z.begin(), z.end(),
            [](const ScoreIdx& a, const ScoreIdx& b) { return a.score > b.score; });

  std::vector<bool> used(na, false);
  std::vector<bool> input_used(hits.size(), false);

  for (std::size_t k = 0; k < z.size() && result.chains.size() < config_.max_chains; ++k) {
    if (used[z[k].idx]) continue;

    std::vector<std::size_t> chain_indices;
    bool shared_prefix = false;
    {
      std::int32_t cur = static_cast<std::int32_t>(z[k].idx);
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

    Chain chain;
    chain.score = static_cast<double>(z[k].score);
    chain.path_id = 0;
    chain.coord_space = CoordSpace::kCanonical;
    chain.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
      auto si = anchors[idx].src_idx;
      const auto& src = hits[si];
      input_used[si] = true;

      ChainedAnchor ca;
      ca.node_id = src.node_id;
      ca.offset = src.offset;
      ca.length = anchors[idx].span;
      ca.read_pos = anchors[idx].query_pos;
      ca.ref_coord = anchors[idx].ref_coord;
      chain.anchors.push_back(ca);
    }

    result.chains.push_back(std::move(chain));
  }

  result.used_inputs = std::move(input_used);
  return result;
}

}  // namespace piru::mapping
