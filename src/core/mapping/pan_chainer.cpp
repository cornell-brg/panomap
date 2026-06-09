/**
 * pan_chainer.cpp
 *
 * PanChainer: 1D-banded cross-path colinear chaining.
 * Sort by canonical 1D coord, band for candidate selection, evaluate
 * colinearity on exact per-path ref coords.
 *
 * Related:
 *  - pan_chainer.hpp
 *  - sort_chainer.cpp (SoA + permutation pattern reused here)
 *  - path_chainer.cpp (per-path scoring logic reused here)
 *
 * SPDX-License-Identifier: MIT
 */

#include "core/mapping/pan_chainer.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "core/util/logging.hpp"

namespace panomap::mapping {

namespace {

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

/* Per-anchor path coordinate entry. Packed for cache friendliness. */
struct PathCoord {
  std::uint32_t path_id;
  std::int64_t ref_coord;
};

/* CSR-like storage for per-anchor path coords + path membership bitmap.
 * anchor i's path coords are in path_coords[offsets[i] .. offsets[i+1]).
 * Sorted by path_id within each anchor for fast intersection.
 * Bitmap: flat array, words_per_anchor uint64 words per anchor.
 * shares_path(i, j) does a quick AND check before merge-intersect. */
struct AnchorPathIndex {
  std::vector<PathCoord> coords;       // flat array of all (path_id, ref_coord)
  std::vector<std::uint32_t> offsets;  // offsets[i] = start of anchor i's coords
  std::vector<std::uint64_t> bitmap;   // path membership bitmap, contiguous SoA
  std::size_t words_per_anchor{0};     // ceil(num_paths / 64)

  /* Quick shared-path check via bitmap AND. */
  inline bool shares_path(std::size_t i, std::size_t j) const {
    auto bi = i * words_per_anchor;
    auto bj = j * words_per_anchor;
    for (std::size_t w = 0; w < words_per_anchor; ++w) {
      if (bitmap[bi + w] & bitmap[bj + w]) return true;
    }
    return false;
  }

  void build(const std::vector<NodeAnchor>& hits,
             const std::vector<std::vector<index::LinearCoordinate>>& lin_coords,
             const std::uint32_t* order, std::size_t na) {
    /* 1. Find max path_id for bitmap sizing */
    std::uint32_t max_pid = 0;
    for (std::size_t k = 0; k < na; ++k) {
      auto nid = hits[order[k]].node_id;
      if (nid < lin_coords.size()) {
        for (const auto& lc : lin_coords[nid])
          if (static_cast<std::uint32_t>(lc.path_id) > max_pid)
            max_pid = static_cast<std::uint32_t>(lc.path_id);
      }
    }
    words_per_anchor = (max_pid / 64) + 1;
    bitmap.assign(na * words_per_anchor, 0);

    /* 2. Build CSR coords + bitmap */
    std::size_t total = 0;
    for (std::size_t k = 0; k < na; ++k) {
      auto nid = hits[order[k]].node_id;
      if (nid < lin_coords.size()) total += lin_coords[nid].size();
    }

    coords.resize(total);
    offsets.resize(na + 1);

    std::size_t pos = 0;
    for (std::size_t k = 0; k < na; ++k) {
      offsets[k] = static_cast<std::uint32_t>(pos);
      auto nid = hits[order[k]].node_id;
      auto ofs = hits[order[k]].offset;
      if (nid < lin_coords.size()) {
        for (const auto& lc : lin_coords[nid]) {
          auto pid = static_cast<std::uint32_t>(lc.path_id);
          coords[pos].path_id = pid;
          coords[pos].ref_coord = lc.ref_coord + static_cast<std::int64_t>(ofs);
          ++pos;
          // Set bitmap bit
          bitmap[k * words_per_anchor + pid / 64] |= (1ULL << (pid % 64));
        }
        std::sort(coords.data() + offsets[k], coords.data() + pos,
                  [](const PathCoord& a, const PathCoord& b) { return a.path_id < b.path_id; });
      }
    }
    offsets[na] = static_cast<std::uint32_t>(pos);
  }
};

/* SoA buffers for canonical 1D coords, query pos, and DP state. */
struct DPBuffers {
  std::vector<double> canon_coord;  // canonical 1D position (for sort/band)
  std::vector<std::uint16_t> query_pos;
  std::vector<std::uint16_t> span;
  std::vector<std::uint32_t> src_idx;  // back-pointer to original hits[]

  std::vector<std::int32_t> q_span;  // pore_k + span - 1

  std::vector<std::int32_t> f, v, pred;
  std::vector<std::uint32_t> chain_path;  // path_id of winning transition ending here
  std::vector<std::uint32_t> order;       // sort permutation

  void resize(std::size_t n) {
    canon_coord.resize(n);
    query_pos.resize(n);
    span.resize(n);
    src_idx.resize(n);
    q_span.resize(n);
    f.resize(n);
    v.resize(n);
    chain_path.resize(n);
    pred.resize(n);
    order.resize(n);
  }
};

void apply_permutation(DPBuffers& dp, std::size_t n) {
  std::vector<double> tmp_f64(n);
  std::vector<std::uint32_t> tmp_u32(n);
  std::vector<std::uint16_t> tmp_u16(n);

  for (std::size_t i = 0; i < n; ++i) tmp_f64[i] = dp.canon_coord[dp.order[i]];
  std::copy(tmp_f64.begin(), tmp_f64.begin() + n, dp.canon_coord.begin());

  for (std::size_t i = 0; i < n; ++i) tmp_u32[i] = dp.src_idx[dp.order[i]];
  std::copy(tmp_u32.begin(), tmp_u32.begin() + n, dp.src_idx.begin());

  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.query_pos[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.query_pos.begin());

  for (std::size_t i = 0; i < n; ++i) tmp_u16[i] = dp.span[dp.order[i]];
  std::copy(tmp_u16.begin(), tmp_u16.begin() + n, dp.span.begin());
}

/* Score transition j->i on a specific shared path.
 * Returns net score (reward - penalty), or INT32_MIN if invalid.
 * reward = min(q_span_j, dg), penalty = gap cost. Same as minimap2/PathChainer. */
inline std::int32_t score_transition(std::int64_t ref_j, std::int64_t ref_i, std::int32_t dq,
                                     std::int32_t q_span_j, std::int32_t max_dist_ref,
                                     std::int32_t bw, float chn_pen_gap, float chn_pen_skip) {
  auto dr = static_cast<std::int32_t>(ref_i - ref_j);
  if (dr <= 0 || dr > max_dist_ref) return std::numeric_limits<std::int32_t>::min();

  auto dd = dr > dq ? dr - dq : dq - dr;
  if (dd > bw) return std::numeric_limits<std::int32_t>::min();

  auto dg = dr < dq ? dr : dq;
  std::int32_t sc = q_span_j < dg ? q_span_j : dg;
  if (dd || dg > q_span_j) {
    float lin_pen = chn_pen_gap * static_cast<float>(dd) + chn_pen_skip * static_cast<float>(dg);
    float log_pen = dd >= 1 ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;
    sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
  }
  return sc;
}

}  // namespace

PanChainer::PanChainer(
    PanChainerConfig config, const std::vector<float>& node_1d_coords,
    const std::vector<std::vector<index::LinearCoordinate>>& linearization_coords,
    const std::vector<std::size_t>& path_lengths)
    : config_(std::move(config)),
      node_1d_coords_(node_1d_coords),
      coords_(linearization_coords),
      path_lengths_(path_lengths) {}

ChainResult PanChainer::chain(const std::vector<NodeAnchor>& hits) const {
  ChainResult result;
  if (hits.empty()) return result;

  /* 1. Convert NodeAnchors to SoA: canonical 1D coord + query pos */

  DPBuffers dp;
  dp.resize(hits.size());

  std::size_t na = 0;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (h.node_id >= node_1d_coords_.size()) continue;
    dp.canon_coord[na] = node_1d_coords_[h.node_id] + static_cast<double>(h.offset);
    dp.query_pos[na] = static_cast<std::uint16_t>(h.read_pos);
    dp.span[na] = h.span;
    dp.src_idx[na] = static_cast<std::uint32_t>(i);
    dp.order[na] = static_cast<std::uint32_t>(na);
    ++na;
  }

  result.expanded_anchor_count = na;
  if (na < config_.min_chain_anchors) return result;

  /* 2. Sort by canonical 1D coord */

  std::sort(dp.order.begin(), dp.order.begin() + na, [&dp](std::uint32_t a, std::uint32_t b) {
    if (dp.canon_coord[a] != dp.canon_coord[b]) return dp.canon_coord[a] < dp.canon_coord[b];
    return dp.query_pos[a] < dp.query_pos[b];
  });
  apply_permutation(dp, na);

  /* 3. Build per-anchor path coord index (CSR-like, sorted by path_id) */

  AnchorPathIndex path_idx;
  path_idx.build(hits, coords_, dp.src_idx.data(), na);

  /* 4. Precompute q_span */

  const bool has_pore_k = config_.pore_k > 0;
  const auto pore_k_i32 = static_cast<std::int32_t>(config_.pore_k);
  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t raw = static_cast<std::int32_t>(dp.span[i]);
    dp.q_span[i] = has_pore_k ? pore_k_i32 + raw - 1 : raw;
  }

  /* 5. Initialize DP */

  std::memset(dp.f.data(), 0, na * sizeof(std::int32_t));
  std::memset(dp.v.data(), 0, na * sizeof(std::int32_t));
  std::fill(dp.pred.begin(), dp.pred.begin() + na, -1);
  std::fill(dp.chain_path.begin(), dp.chain_path.begin() + na, UINT32_MAX);

  /* 6. DP forward pass: 1D band + shared-path scoring */

  const double band_1d = static_cast<double>(config_.band_1d);
  const auto max_dist_ref = static_cast<std::int32_t>(config_.max_dist_ref);
  const auto max_dist_query = static_cast<std::int32_t>(config_.max_dist_query);
  const auto bw = static_cast<std::int32_t>(config_.bw);
  const float chn_pen_gap = config_.chn_pen_gap;
  const float chn_pen_skip = config_.chn_pen_skip;
  const auto chn_pen_switch = static_cast<std::int32_t>(config_.chn_pen_switch);

  const auto* cc = dp.canon_coord.data();
  const auto* qp = dp.query_pos.data();
  const auto* qs = dp.q_span.data();
  auto* f = dp.f.data();
  auto* v = dp.v.data();
  auto* pred = dp.pred.data();
  auto* chain_path = dp.chain_path.data();

  std::size_t st = 0;

  for (std::size_t i = 0; i < na; ++i) {
    std::int32_t max_f = qs[i];
    std::int32_t max_j = -1;

    // Advance 1D band window
    while (st < i && cc[i] - cc[st] > band_1d) ++st;

    // Scan backwards within 1D band
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

      // Quick bitmap check: skip if no shared path possible
      if (!path_idx.shares_path(i, j)) continue;

      /* Shared-path intersection + scoring via merge of sorted path_id lists */
      std::int32_t best_net = std::numeric_limits<std::int32_t>::min();
      std::uint32_t best_pid = UINT32_MAX;
      {
        auto i_begin = path_idx.offsets[i], i_end = path_idx.offsets[i + 1];
        auto j_begin = path_idx.offsets[j], j_end = path_idx.offsets[j + 1];
        auto ip = i_begin, jp = j_begin;

        while (ip < i_end && jp < j_end) {
          auto pid_i = path_idx.coords[ip].path_id;
          auto pid_j = path_idx.coords[jp].path_id;
          if (pid_i < pid_j) {
            ++ip;
            continue;
          }
          if (pid_j < pid_i) {
            ++jp;
            continue;
          }
          // Shared path found -- score transition (returns net = reward - penalty)
          auto net = score_transition(path_idx.coords[jp].ref_coord, path_idx.coords[ip].ref_coord,
                                      dq, qs[j], max_dist_ref, bw, chn_pen_gap, chn_pen_skip);
          // Subtract switch penalty if chain at j was on a different path
          if (net > std::numeric_limits<std::int32_t>::min() && chain_path[j] != UINT32_MAX &&
              chain_path[j] != pid_i) {
            net -= chn_pen_switch;
          }
          if (net > best_net) {
            best_net = net;
            best_pid = pid_i;
          }
          ++ip;
          ++jp;
        }
      }

      if (best_net == std::numeric_limits<std::int32_t>::min()) {
        continue;  // no shared path -- don't count toward skip budget
      }

      std::int32_t sc = best_net + f[j];
      if (sc > max_f) {
        max_f = sc;
        max_j = static_cast<std::int32_t>(j);
        chain_path[i] = best_pid;
        if (num_skipped > 0) --num_skipped;
      } else {
        ++num_skipped;
      }
    }

    f[i] = max_f;
    pred[i] = max_j;
    // chain_path[i] already set during inner loop (or stays UINT32_MAX if new chain)
    v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
  }

  /* 7. Multi-chain extraction */

  std::vector<bool> used(na, false);
  std::vector<bool> input_used(hits.size(), false);

  while (result.chains.size() < config_.max_chains) {
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

    // Backtrack
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

    /* Determine best path for the chain.
     * Walk chain, find the path that appears most consistently. */
    Chain chain;
    chain.score = static_cast<double>(best_score);
    chain.coord_space = CoordSpace::kPath;
    chain.anchors.reserve(chain_indices.size());

    // Find dominant path: for each consecutive pair, count which path_ids work
    std::vector<std::uint32_t> path_votes;
    for (std::size_t k = 0; k + 1 < chain_indices.size(); ++k) {
      auto ci = chain_indices[k], ci_next = chain_indices[k + 1];
      auto a_begin = path_idx.offsets[ci], a_end = path_idx.offsets[ci + 1];
      auto b_begin = path_idx.offsets[ci_next], b_end = path_idx.offsets[ci_next + 1];
      auto ap = a_begin, bp = b_begin;
      while (ap < a_end && bp < b_end) {
        auto pa = path_idx.coords[ap].path_id;
        auto pb = path_idx.coords[bp].path_id;
        if (pa < pb) {
          ++ap;
          continue;
        }
        if (pb < pa) {
          ++bp;
          continue;
        }
        path_votes.push_back(pa);
        ++ap;
        ++bp;
      }
    }

    // Pick most frequent path_id
    std::uint32_t best_path = 0;
    if (!path_votes.empty()) {
      std::sort(path_votes.begin(), path_votes.end());
      std::uint32_t cur_pid = path_votes[0], cur_count = 1, best_count = 0;
      for (std::size_t k = 1; k < path_votes.size(); ++k) {
        if (path_votes[k] == cur_pid) {
          ++cur_count;
        } else {
          if (cur_count > best_count) {
            best_count = cur_count;
            best_path = cur_pid;
          }
          cur_pid = path_votes[k];
          cur_count = 1;
        }
      }
      if (cur_count > best_count) best_path = cur_pid;
    }

    chain.path_id = best_path;

    // Build chain anchors using the dominant path's ref coords
    for (std::size_t idx : chain_indices) {
      auto si = dp.src_idx[idx];
      const auto& src = hits[si];
      input_used[si] = true;

      ChainedAnchor ca;
      ca.node_id = src.node_id;
      ca.offset = src.offset;
      ca.length = dp.span[idx];
      ca.read_pos = dp.query_pos[idx];

      // Look up ref_coord on the dominant path
      ca.ref_coord = static_cast<std::int64_t>(dp.canon_coord[idx]);
      auto a_begin = path_idx.offsets[idx], a_end = path_idx.offsets[idx + 1];
      for (auto p = a_begin; p < a_end; ++p) {
        if (path_idx.coords[p].path_id == best_path) {
          ca.ref_coord = path_idx.coords[p].ref_coord;
          break;
        }
      }

      chain.anchors.push_back(ca);
    }

    result.chains.push_back(std::move(chain));
  }

  result.used_inputs = std::move(input_used);
  return result;
}

/* PanChainerConfig CLI integration */

PanChainerConfig PanChainerConfig::from_parsed(const cli::Parsed& parsed) {
  PanChainerConfig cfg;
  if (parsed.values.count("chain-band-1d"))
    cfg.band_1d = std::stoull(parsed.values.at("chain-band-1d"));
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
  if (parsed.values.count("chain-pen-switch"))
    cfg.chn_pen_switch = std::stof(parsed.values.at("chain-pen-switch"));
  return cfg;
}

}  // namespace panomap::mapping
