// SPDX-License-Identifier: MIT

#include "io/results/paf_writer.hpp"

#include <cstdio>
#include <sstream>

#include "util/logging.hpp"

namespace piru::io {

PafWriter::PafWriter(const std::string& path) : out_(path, std::ios::out | std::ios::trunc) {
  if (!out_) {
    LOG_ERROR("Failed to open PAF output: " + path);
  }
}

PafWriter::~PafWriter() { out_.flush(); }

void PafWriter::write(const mapping::ReadMapResult& result, const std::string& read_id,
                      std::size_t read_length) {
  if (!out_) return;

  /* Unmapped */
  if (!result.mapped()) {
    std::stringstream ss;
    if (!result.mappings.empty()) {
      const auto& best = result.mappings[0];
      std::uint32_t min_q = best.anchors.front().read_pos;
      std::uint32_t max_q = best.anchors.front().read_pos;
      std::int64_t min_r = best.anchors.front().ref_coord;
      std::int64_t max_r = best.anchors.front().ref_coord;
      for (const auto& a : best.anchors) {
        min_q = std::min(min_q, a.read_pos);
        max_q = std::max(max_q, a.read_pos + a.length);
        min_r = std::min(min_r, a.ref_coord);
        max_r = std::max(max_r, a.ref_coord + static_cast<std::int64_t>(a.length));
      }
      ss << read_id << '\t' << read_length << '\t' << min_q << '\t' << max_q
         << "\t*\t*\t0\t*\t*\t0\t0\t" << best.mapq;
      ss << "\ttp:A:U";
      ss << "\tcs:i:" << static_cast<int>(best.chain_score);
      ss << "\tan:i:" << best.anchors.size();
      auto qspan = max_q - min_q;
      if (qspan > 0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "se:f:%.3f", best.chain_score / static_cast<double>(qspan));
        ss << '\t' << buf;
      }
      auto rspan = max_r - min_r;
      if (rspan > 0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "ad:f:%.4f",
                      static_cast<double>(best.anchors.size()) / static_cast<double>(rspan));
        ss << '\t' << buf;
      }
    } else {
      ss << read_id << '\t' << read_length << "\t*\t*\t*\t*\t0\t*\t*\t0\t0\t0\ttp:A:U";
    }
    {
      char dt_buf[32];
      std::snprintf(dt_buf, sizeof(dt_buf), "dt:f:%.6f", result.processing_time_sec);
      ss << '\t' << dt_buf;
    }
    out_ << ss.str() << '\n';
    return;
  }

  /* Mapped: primary */
  const auto& mapping = result.mappings[0];

  std::size_t min_query = mapping.anchors.front().read_pos;
  std::size_t max_query = mapping.anchors.front().read_pos;
  for (const auto& a : mapping.anchors) {
    min_query = std::min(min_query, static_cast<std::size_t>(a.read_pos));
    max_query = std::max(max_query, static_cast<std::size_t>(a.read_pos) + a.length);
  }

  std::int64_t min_ref = mapping.anchors.front().ref_coord;
  std::int64_t max_ref = mapping.anchors.front().ref_coord;
  for (const auto& a : mapping.anchors) {
    min_ref = std::min(min_ref, a.ref_coord);
    max_ref = std::max(max_ref, a.ref_coord + static_cast<std::int64_t>(a.length));
  }

  std::uint64_t total_anchor_len = 0;
  for (const auto& a : mapping.anchors) total_anchor_len += a.length;

  auto qspan = max_query - min_query;

  std::stringstream ss;
  ss << read_id << '\t' << read_length << '\t' << min_query << '\t' << max_query << "\t+\t*\t0\t"
     << min_ref << '\t' << max_ref << '\t' << total_anchor_len << '\t' << (max_ref - min_ref)
     << '\t' << mapping.mapq;
  ss << "\ttp:A:P";
  ss << "\tcs:i:" << static_cast<int>(mapping.chain_score);
  ss << "\tan:i:" << mapping.anchors.size();
  if (qspan > 0) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "se:f:%.3f", mapping.chain_score / static_cast<double>(qspan));
    ss << '\t' << buf;
  }

  out_ << ss.str() << '\n';
}

}  // namespace piru::io
