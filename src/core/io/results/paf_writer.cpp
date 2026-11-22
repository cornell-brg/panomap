// SPDX-License-Identifier: MIT

#include "core/io/results/paf_writer.hpp"

#include <algorithm>
#include <sstream>
#include <string_view>

#include "core/io/results/result_tags.hpp"
#include "core/util/logging.hpp"

namespace panomap::io {

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
  std::stringstream ss;

  /* Unmapped: no target; keep the query span and the stats of the best chain, if any. */
  if (!result.mapped()) {
    if (!result.mappings.empty()) {
      const auto& best = result.mappings[0];
      const auto s = tags::chain_span(best);
      ss << read_id << '\t' << read_length << '\t' << s.q_lo << '\t' << s.q_hi
         << "\t*\t*\t0\t*\t*\t0\t0\t" << best.mapq;
      tags::write_chain(ss, 'U', best, s);
    } else {
      ss << read_id << '\t' << read_length << "\t*\t*\t*\t*\t0\t*\t*\t0\t0\t0\ttp:A:U";
    }
    tags::write_read(ss, result);
    out_ << ss.str() << '\n';
    return;
  }

  /* Mapped: primary only. The target is the chosen path. A path named "<name>_reverse"
   * is the reverse strand of <name>, with coordinates flipped. */
  const auto& mapping = result.mappings[0];
  const auto s = tags::chain_span(mapping);
  std::string target = "*";
  std::int64_t target_len = 0;
  char strand = '+';
  std::int64_t min_ref = s.r_lo;
  std::int64_t max_ref = s.r_hi;
  const bool on_path = mapping.coord_space != mapping::CoordSpace::kCanonical &&
                       mapping.path_id < graph_.pathCount();
  if (on_path) {
    const auto pid = static_cast<std::uint32_t>(mapping.path_id);
    target = std::string(graph_.pathName(pid));
    target_len = static_cast<std::int64_t>(graph_.pathLength(pid));
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
  const std::int64_t block = max_ref - min_ref;

  ss << read_id << '\t' << read_length << '\t' << s.q_lo << '\t' << s.q_hi << '\t' << strand
     << '\t' << target << '\t' << target_len << '\t' << min_ref << '\t' << max_ref << '\t'
     << tags::matches(mapping, block) << '\t' << block << '\t' << mapping.mapq;
  tags::write_chain(ss, 'P', mapping, s);
  tags::write_read(ss, result);
  out_ << ss.str() << '\n';
}

}  // namespace panomap::io
