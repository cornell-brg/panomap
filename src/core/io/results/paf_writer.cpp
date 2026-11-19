// SPDX-License-Identifier: MIT

#include "core/io/results/paf_writer.hpp"

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <string_view>

#include "core/util/logging.hpp"

namespace panomap::io {

namespace {
// Timing tags in the form common PAF scoring scripts read: mt:f = map time in ms,
// ci:i = number of chunks processed.
void write_timing_tags(std::ostream& ss, const mapping::ReadMapResult& result) {
  char mt_buf[32];
  std::snprintf(mt_buf, sizeof(mt_buf), "mt:f:%.3f", result.processing_time_sec * 1000.0);
  ss << '\t' << mt_buf << "\tci:i:" << result.chunks_processed;
}
}  // namespace

PafWriter::PafWriter(const std::string& path, const index::FlatGraph& graph)
    : out_(path, std::ios::out | std::ios::trunc), graph_(graph) {
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
    write_timing_tags(ss, result);
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

  /* Target name, length, span and strand: same rules as the GAF writer. A path named
   * "<name>_reverse" is the reverse strand of <name>, with coordinates flipped. */
  std::string target = "*";
  std::int64_t target_len = 0;
  char strand = '+';
  std::int64_t min_ref = mapping.anchors.front().ref_coord;
  std::int64_t max_ref = mapping.anchors.front().ref_coord;
  const bool on_path = mapping.coord_space != mapping::CoordSpace::kCanonical &&
                       mapping.path_id < graph_.pathCount();
  if (on_path) {
    const auto pid = static_cast<std::uint32_t>(mapping.path_id);
    target = std::string(graph_.pathName(pid));
    target_len = static_cast<std::int64_t>(graph_.pathLength(pid));
  }
  for (const auto& a : mapping.anchors) {
    const auto end = a.ref_coord + static_cast<std::int64_t>(a.length);
    min_ref = std::min(min_ref, a.ref_coord);
    max_ref = std::max(max_ref, end);
  }
  if (on_path) {
    min_ref = std::max(std::int64_t{0}, min_ref);
    if (target_len > 0) max_ref = std::min(max_ref, target_len);
    constexpr std::string_view kRevSuffix = "_reverse";
    if (target.size() > kRevSuffix.size() &&
        std::string_view(target).substr(target.size() - kRevSuffix.size()) == kRevSuffix) {
      strand = '-';
      target.resize(target.size() - kRevSuffix.size());
      const auto flipped_start = target_len - max_ref;
      max_ref = target_len - min_ref;
      min_ref = flipped_start;
    }
  }

  std::uint64_t total_anchor_len = 0;
  for (const auto& a : mapping.anchors) total_anchor_len += a.length;

  auto qspan = max_query - min_query;

  std::stringstream ss;
  ss << read_id << '\t' << read_length << '\t' << min_query << '\t' << max_query << '\t' << strand
     << '\t' << target << '\t' << target_len << '\t' << min_ref << '\t' << max_ref << '\t' << total_anchor_len << '\t' << (max_ref - min_ref)
     << '\t' << mapping.mapq;
  ss << "\ttp:A:P";
  ss << "\tcs:i:" << static_cast<int>(mapping.chain_score);
  ss << "\tan:i:" << mapping.anchors.size();
  if (qspan > 0) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "se:f:%.3f", mapping.chain_score / static_cast<double>(qspan));
    ss << '\t' << buf;
  }
  {
    char ws_buf[32];
    std::snprintf(ws_buf, sizeof(ws_buf), "ws:f:%.3f", result.standout);
    ss << '\t' << ws_buf;
    ss << "\tnc:i:" << result.mappings.size();
    ss << "\tdf:A:" << static_cast<char>(result.decision_path);
    ss << "\tck:i:" << result.chunks_processed;
    char dt_buf[32];
    std::snprintf(dt_buf, sizeof(dt_buf), "dt:f:%.6f", result.processing_time_sec);
    ss << '\t' << dt_buf;
  }
  write_timing_tags(ss, result);

  out_ << ss.str() << '\n';
}

}  // namespace panomap::io
