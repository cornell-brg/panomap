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

#include "core/mapping/sort_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "core/util/logging.hpp"

namespace panomap::mapping {

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
 * loop touches only the arrays it needs (ref_coord, query_pos, q_span).
 *
 * Type choices:
 *  ref_coord: uint64 -- 1D canonical coords, always positive. Can exceed 4.2B
 *             for whole-genome pangenomes with many disconnected components
 *             (gaps between components accumulate on the 1D axis).
 *  query_pos: uint16 -- event positions in read, max ~4800 at 10 chunks
 *  span:      uint16 -- anchor coverage length (from NodeAnchor)
 *  q_span:    int32  -- pore_k + span - 1, can exceed uint16 after merge */
struct DPBuffers {
  /* Anchor data */
  std::vector<double> ref_coord;
  std::vector<std::uint16_t> query_pos;
  std::vector<std::uint16_t> span;
  std::vector<std::uint32_t> src_idx;
  std::vector<std::uint32_t> cc_id;  // connected component of the anchor's node

  /* Precomputed per-anchor */
  std::vector<std::int32_t> q_span;

  /* DP state */
  std::vector<std::int32_t> f, v, pred, t;
  std::vector<float> ratio;  // dr/dq of last transition in best chain ending here

  /* Sort permutation */
  std::vector<std::uint32_t> order;

  void resize(std::size_t n) {
    ref_coord.resize(n);
    query_pos.resize(n);
    span.resize(n);
    src_idx.resize(n);
    cc_id.resize(n);
    q_span.resize(n);
    f.resize(n);
    v.resize(n);
    pred.resize(n);
    t.resize(n);
    ratio.resize(n);
    order.resize(n);
  }
};

// Apply sort permutation in-place to all SoA arrays.
void apply_permutation(DPBuffers& dp, std::size_t n) {
  std::vector<double> tmp_f64(n);
  std::vector<std::uint32_t> tmp_u32(n);
  std::vector<std::uint16_t> tmp_u16(n);

  /* ref_coord (double) */
  for (std::size_t i = 0; i < n; ++i) tmp_f64[i] = dp.ref_coord[dp.order[i]];
  std::copy(tmp_f64.begin(), tmp_f64.begin() + n, dp.ref_coord.begin());

  /* src_idx (uint32) */
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.src_idx[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.src_idx.begin());

  /* query_pos (uint16) */
  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.query_pos[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.query_pos.begin());

  /* span (uint16) */
  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.span[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.span.begin());

  /* cc_id (uint32) */
  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.cc_id[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.cc_id.begin());
}

}  // namespace

SortChainer::SortChainer(SortChainerConfig config, const std::vector<float>& node_1d_coords,
                         std::vector<std::uint32_t> node_bp_lens,
                         const std::vector<std::uint32_t>& component_ids)
    : config_(std::move(config)),
      node_1d_coords_(node_1d_coords),
      node_bp_lens_(std::move(node_bp_lens)),
      component_ids_(component_ids) {}

ChainResult SortChainer::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult fwd = chain_pass(hits, /*reverse_dir=*/false);
  if (!config_.bidirectional) return fwd;

  ChainResult bwd = chain_pass(hits, /*reverse_dir=*/true);

  /* Merge: combine chains from both passes, sort by score desc, take top-K. */
  ChainResult merged;
  merged.expanded_anchor_count = std::max(fwd.expanded_anchor_count, bwd.expanded_anchor_count);
  merged.chains.reserve(fwd.chains.size() + bwd.chains.size());
  for (auto& c : fwd.chains) merged.chains.push_back(std::move(c));
  for (auto& c : bwd.chains) merged.chains.push_back(std::move(c));
  std::sort(merged.chains.begin(), merged.chains.end(),
            [](const Chain& a, const Chain& b) { return a.score > b.score; });
  if (merged.chains.size() > config_.max_chains) {
    merged.chains.resize(config_.max_chains);
  }
  /* used_inputs: logical OR of both passes (anchor used by either pass). */
  merged.used_inputs.assign(hits.size(), false);
  for (std::size_t i = 0; i < hits.size(); ++i) {
    bool a = i < fwd.used_inputs.size() && fwd.used_inputs[i];
    bool b = i < bwd.used_inputs.size() && bwd.used_inputs[i];
    merged.used_inputs[i] = a || b;
  }
  return merged;
}

ChainResult SortChainer::chain_pass(const std::vector<NodeAnchor>& hits, bool reverse_dir) const {
  ChainResult result;
  if (hits.empty()) return result;

  /* 1. Convert NodeAnchors to SoA arrays.
   * For reverse_dir, negate ref_coord so the same DP sort+scan logic finds
   * chains where original-1D-coord decreases as query_pos increases. */

  DPBuffers dp;
  dp.resize(hits.size());

  const bool has_cc = !component_ids_.empty();
  std::size_t na = 0;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (h.node_id >= node_1d_coords_.size()) continue;
    // Scale offset by node's 1D-to-bp ratio (forward=even, reverse/end=odd)
    double node_start = node_1d_coords_[h.node_id];
    double node_end = node_1d_coords_[h.node_id | 1];
    double node_1d_span = node_end - node_start;
    double ref;
    if (!node_bp_lens_.empty() && h.node_id < node_bp_lens_.size() &&
        node_bp_lens_[h.node_id] > 0 && node_1d_span > 0) {
      ref = node_start + static_cast<double>(h.offset) * (node_1d_span / node_bp_lens_[h.node_id]);
    } else {
      ref = node_start + static_cast<double>(h.offset);
    }
    if (reverse_dir) ref = -ref;
    dp.ref_coord[na] = ref;
    dp.query_pos[na] = static_cast<std::uint16_t>(h.read_pos);
    dp.span[na] = h.span;
    dp.src_idx[na] = static_cast<std::uint32_t>(i);
    // Look up CC by node id. If component_ids empty, treat as single CC=0.
    dp.cc_id[na] = has_cc && h.node_id < component_ids_.size() ? component_ids_[h.node_id] : 0u;
    dp.order[na] = static_cast<std::uint32_t>(na);
    ++na;
  }

  result.expanded_anchor_count = na;
  if (na < config_.min_chain_anchors) return result;

  /* 2. Sort by (cc_id, ref_coord, query_pos) via index permutation.
   *    CC-first keeps each CC contiguous so the outer DP loop can reset its
   *    window at CC boundaries (preventing cross-CC chaining). */

  std::sort(dp.order.begin(), dp.order.begin() + na, [&dp](std::uint32_t a, std::uint32_t b) {
    if (dp.cc_id[a] != dp.cc_id[b]) return dp.cc_id[a] < dp.cc_id[b];
    if (dp.ref_coord[a] != dp.ref_coord[b]) return dp.ref_coord[a] < dp.ref_coord[b];
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

  std::memset(dp.f.data(), 0, na * sizeof(std::int32_t));
  std::memset(dp.v.data(), 0, na * sizeof(std::int32_t));
  std::memset(dp.t.data(), 0, na * sizeof(std::int32_t));
  std::fill(dp.pred.begin(), dp.pred.begin() + na, -1);
  std::fill(dp.ratio.begin(), dp.ratio.begin() + na, 0.0f);

  /* 5. DP forward pass */

  const auto max_dist_ref = static_cast<double>(config_.max_dist_ref);
  const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
  const auto bw = static_cast<std::int32_t>(config_.bw);
  const float chn_pen_gap = config_.chn_pen_gap;
  const float chn_pen_skip = config_.chn_pen_skip;
  const float dd_tol_frac = config_.dd_tolerance_frac;
  const float chn_pen_ratio = config_.chn_pen_ratio;

  // Raw pointers for inner loop (keeps them in registers)
  const auto* rc = dp.ref_coord.data();
  const auto* qp = dp.query_pos.data();
  const auto* qs = dp.q_span.data();
  const auto* cc = dp.cc_id.data();
  auto* f = dp.f.data();
  auto* v = dp.v.data();
  auto* pred = dp.pred.data();
  auto* t = dp.t.data();
  auto* ratio = dp.ratio.data();

  std::size_t st = 0;
  std::uint32_t prev_cc = ~0u;  // sentinel: no CC seen yet
  std::int64_t max_ii = -1;

  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t max_f = qs[i];
    std::int32_t max_j = -1;

    // CC boundary: anchors are sorted by (cc_id, ref_coord), so when we cross
    // into a new CC, all earlier anchors are in different CCs. Reset st to i
    // so the predecessor scan only considers same-CC anchors.
    if (cc[i] != prev_cc) {
      st = i;
      prev_cc = cc[i];
    }
    // Advance start pointer past anchors too far back in ref (same CC guaranteed)
    while (st < i && rc[i] > rc[st] + max_dist_ref) ++st;

    /* Scan backwards for predecessor candidates */
    std::size_t num_skipped = 0;
    std::size_t num_iter = 0;
    const std::size_t max_iter = config_.max_iterations;
    for (std::size_t j = i; j > st && num_skipped < config_.max_skip;) {
      --j;
      if (max_iter > 0 && ++num_iter > max_iter) break;

      auto dq = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[j]);
      if (dq <= 0 || dq > max_dist_query) {
        ++num_skipped;
        continue;
      }

      double dr_d = rc[i] - rc[j];
      if (dr_d <= 0.0 || dr_d > max_dist_ref) {
        ++num_skipped;
        continue;
      }
      auto dr = static_cast<std::int32_t>(dr_d + 0.5);  // round to nearest

      auto dd_raw = dr > dq ? dr - dq : dq - dr;
      if (dd_raw > bw) {
        ++num_skipped;
        continue;
      }

      /* Compute score for j -> i transition */
      auto dg = dr < dq ? dr : dq;
      // Dead zone: dd within tolerance_frac * dg is penalty-free
      auto dd_tol = static_cast<std::int32_t>(dd_tol_frac * static_cast<float>(dg));
      auto dd = dd_raw > dd_tol ? dd_raw - dd_tol : 0;
      std::int32_t sc = qs[j] < dg ? qs[j] : dg;
      if (dd || dg > qs[j]) {
        float lin_pen =
            chn_pen_gap * static_cast<float>(dd) + chn_pen_skip * static_cast<float>(dg);
        float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
        sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
      }
      // Ratio consistency: penalize if dr/dq deviates from chain's running ratio
      if (chn_pen_ratio > 0.0f && pred[j] >= 0) {
        float ratio_ji = static_cast<float>(dr) / static_cast<float>(dq);
        float ratio_dev = ratio_ji > ratio[j] ? ratio_ji - ratio[j] : ratio[j] - ratio_ji;
        sc -= static_cast<std::int32_t>(chn_pen_ratio * ratio_dev * static_cast<float>(dg));
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
        if (f[j] > best) {
          best = f[j];
          max_ii = static_cast<std::int64_t>(j);
        }
      }
    }

    if (max_ii >= 0 && static_cast<std::size_t>(max_ii) >= st) {
      auto dq_m = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[max_ii]);
      double dr_m_d = rc[i] - rc[max_ii];
      std::int32_t tmp = std::numeric_limits<std::int32_t>::min();
      if (dq_m > 0 && dq_m <= max_dist_query && dr_m_d > 0.0 && dr_m_d <= max_dist_ref) {
        auto dr_m = static_cast<std::int32_t>(dr_m_d + 0.5);
        auto dd_m_raw = dr_m > dq_m ? dr_m - dq_m : dq_m - dr_m;
        if (dd_m_raw <= bw) {
          auto dg_m = dr_m < dq_m ? dr_m : dq_m;
          auto dd_m_tol = static_cast<std::int32_t>(dd_tol_frac * static_cast<float>(dg_m));
          auto dd_m = dd_m_raw > dd_m_tol ? dd_m_raw - dd_m_tol : 0;
          tmp = qs[max_ii] < dg_m ? qs[max_ii] : dg_m;
          if (dd_m || dg_m > qs[max_ii]) {
            float lin_pen =
                chn_pen_gap * static_cast<float>(dd_m) + chn_pen_skip * static_cast<float>(dg_m);
            float log_pen = dd_m >= 1 ? mg_log2(static_cast<float>(dd_m + 1)) : 0.0f;
            tmp -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
          }
          if (chn_pen_ratio > 0.0f && pred[max_ii] >= 0) {
            float ratio_mi = static_cast<float>(dr_m) / static_cast<float>(dq_m);
            float ratio_dev =
                ratio_mi > ratio[max_ii] ? ratio_mi - ratio[max_ii] : ratio[max_ii] - ratio_mi;
            tmp -= static_cast<std::int32_t>(chn_pen_ratio * ratio_dev * static_cast<float>(dg_m));
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
    // Store running ratio for this chain
    if (max_j >= 0) {
      auto dq_i = static_cast<std::int32_t>(qp[i]) - static_cast<std::int32_t>(qp[max_j]);
      double dr_i = rc[i] - rc[max_j];
      ratio[i] = (dq_i > 0) ? static_cast<float>(dr_i) / static_cast<float>(dq_i) : 0.0f;
    }
    v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
    if (max_ii < 0 || (rc[i] <= rc[max_ii] + max_dist_ref && f[max_ii] < f[i])) {
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
      if (!used[i] && f[i] >= static_cast<std::int32_t>(config_.min_chain_score) && v[i] <= f[i] &&
          f[i] > best_score) {
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
      // Restore original-sign ref_coord (we negated it for reverse_dir DP).
      ca.ref_coord = reverse_dir ? -dp.ref_coord[idx] : dp.ref_coord[idx];
      chain.anchors.push_back(ca);
    }

    result.chains.push_back(std::move(chain));
  }

  result.used_inputs = std::move(input_used);
  return result;
}

/* SortChainerConfig CLI integration */

SortChainerConfig SortChainerConfig::from_parsed(const cli::Parsed& parsed) {
  SortChainerConfig cfg;
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
  if (parsed.values.count("chain-max-chains"))
    cfg.max_chains = std::stoull(parsed.values.at("chain-max-chains"));
  if (parsed.values.count("chain-max-skip"))
    cfg.max_skip = std::stoull(parsed.values.at("chain-max-skip"));
  if (parsed.values.count("chain-max-iter"))
    cfg.max_iterations = std::stoull(parsed.values.at("chain-max-iter"));
  if (parsed.values.count("chain-dd-tolerance"))
    cfg.dd_tolerance_frac = std::stof(parsed.values.at("chain-dd-tolerance"));
  if (parsed.values.count("chain-pen-ratio"))
    cfg.chn_pen_ratio = std::stof(parsed.values.at("chain-pen-ratio"));
  if (parsed.values.count("chain-bidirectional"))
    cfg.bidirectional = true;
  return cfg;
}

}  // namespace panomap::mapping
