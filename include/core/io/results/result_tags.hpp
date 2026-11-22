/**
 * result_tags.hpp
 *
 * The optional tags that both result writers (PAF and GAF) emit, written in one
 * place so the two formats name and format each fact the same way. The README
 * section "Output tags" lists them.
 *
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <ostream>
#include <utility>
#include <vector>

#include "core/mapping/map_result.hpp"

namespace panomap::io::tags {

// Query span (event positions) and reference span (path or 1D coordinates) of a chain.
struct Span {
  std::int64_t q_lo{0}, q_hi{0}, r_lo{0}, r_hi{0};
};

inline Span chain_span(const mapping::Mapping& m) {
  const auto& a0 = m.anchors.front();
  Span s{a0.read_pos, a0.read_pos, a0.ref_coord, a0.ref_coord};
  for (const auto& a : m.anchors) {
    s.q_lo = std::min<std::int64_t>(s.q_lo, a.read_pos);
    s.q_hi = std::max<std::int64_t>(s.q_hi, std::int64_t{a.read_pos} + a.length);
    s.r_lo = std::min(s.r_lo, a.ref_coord);
    s.r_hi = std::max(s.r_hi, a.ref_coord + static_cast<std::int64_t>(a.length));
  }
  return s;
}

// Reference bases covered by at least one anchor (anchors can overlap), capped at the
// block length, so that matches never exceed the block length.
inline std::int64_t matches(const mapping::Mapping& m, std::int64_t block_len) {
  std::vector<std::pair<std::int64_t, std::int64_t>> iv;
  iv.reserve(m.anchors.size());
  for (const auto& a : m.anchors) iv.emplace_back(a.ref_coord, a.ref_coord + a.length);
  std::sort(iv.begin(), iv.end());
  std::int64_t covered = 0, cur_lo = 0, cur_hi = -1;
  for (const auto& [lo, hi] : iv) {
    if (lo > cur_hi) {
      if (cur_hi > cur_lo) covered += cur_hi - cur_lo;
      cur_lo = lo;
      cur_hi = hi;
    } else {
      cur_hi = std::max(cur_hi, hi);
    }
  }
  if (cur_hi > cur_lo) covered += cur_hi - cur_lo;
  return std::min(covered, std::max<std::int64_t>(block_len, 0));
}

inline void put_f(std::ostream& os, const char* fmt, double v) {
  char buf[48];
  std::snprintf(buf, sizeof(buf), fmt, v);
  os << '\t' << buf;
}

// Chain tags: tp:A (P primary, S secondary, U unmapped), s1:i chain score (RawHash2's
// name; cs is the standard difference-string tag, which vg and minimap2 read),
// an:i anchors, se:f score per query event, ad:f anchors per reference base.
inline void write_chain(std::ostream& os, char type, const mapping::Mapping& m, const Span& s) {
  os << "\ttp:A:" << type << "\ts1:i:" << static_cast<int>(m.chain_score)
     << "\tan:i:" << m.anchors.size();
  if (s.q_hi > s.q_lo) put_f(os, "se:f:%.3f", m.chain_score / static_cast<double>(s.q_hi - s.q_lo));
  if (s.r_hi > s.r_lo) {
    put_f(os, "ad:f:%.4f",
          static_cast<double>(m.anchors.size()) / static_cast<double>(s.r_hi - s.r_lo));
  }
}

// Read tags, once per read: ws:f weighted standout, nc:i chains kept, df:A decision
// path (G gate, S standout, F fallback, U unmapped), mt:f map time in ms, ci:i chunks
// processed. mt and ci are the names RawHash2 and its scripts use.
inline void write_read(std::ostream& os, const mapping::ReadMapResult& r) {
  put_f(os, "ws:f:%.3f", r.standout);
  os << "\tnc:i:" << r.mappings.size() << "\tdf:A:" << static_cast<char>(r.decision_path);
  put_f(os, "mt:f:%.3f", r.processing_time_sec * 1000.0);
  os << "\tci:i:" << r.chunks_processed;
}

}  // namespace panomap::io::tags
