// SPDX-License-Identifier: MIT

#include "io/results/paf_writer.hpp"

#include <sstream>

#include "util/logging.hpp"

namespace piru::io {

PafWriter::PafWriter(const std::string& path) : out_(path, std::ios::out | std::ios::trunc) {
  if (!out_) {
    LOG_ERROR("Failed to open PAF output: " + path);
  }
}

PafWriter::~PafWriter() { out_.flush(); }

void PafWriter::write(const mapping::ReadMapResult& result,
                      const std::string& read_id,
                      std::size_t read_length) {
  if (!out_ || result.mappings.empty()) return;

  // PAF: write primary mapping only (simple for now)
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

  std::stringstream ss;
  ss << read_id << '\t' << read_length << '\t' << min_query << '\t' << max_query
     << "\t+\t*\t0\t" << min_ref << '\t' << max_ref
     << '\t' << total_anchor_len << '\t' << (max_ref - min_ref) << '\t' << mapping.mapq;
  ss << "\tcs:i:" << static_cast<int>(mapping.chain_score);

  out_ << ss.str() << '\n';
}

}  // namespace piru::io
