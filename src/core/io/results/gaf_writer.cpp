/**
 * gaf_writer.cpp
 *
 * Writes mapping results in GAF. For a chain on a path, column 6 is the walk of every
 * path step the chain covers, columns 7 to 9 are the walk length and the chain start
 * and end on the walk, and pn:Z names the path. With anchor_nodes set, column 6 lists
 * only the nodes the anchors fall on, and columns 7 to 9 hold path positions as in PAF.
 * Tags shared with PAF come from result_tags.hpp.
 *
 * SPDX-License-Identifier: MIT
 */

#include "core/io/results/gaf_writer.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <string_view>

#include "core/io/results/result_tags.hpp"
#include "core/util/logging.hpp"

namespace panomap::io {

GafWriter::GafWriter(const std::string& path, const index::FlatGraph& graph, GafWriterConfig config)
    : file_(path, std::ios::out | std::ios::trunc), out_(&file_), graph_(graph), config_(std::move(config)) {
  if (!file_) {
    LOG_ERROR("Failed to open GAF output: " + path);
  }
}

GafWriter::GafWriter(const index::FlatGraph& graph, GafWriterConfig config)
    : out_(&std::cout), graph_(graph), config_(std::move(config)) {}

GafWriter::~GafWriter() { out_->flush(); }

void GafWriter::write(const mapping::ReadMapResult& result, const std::string& read_id,
                      std::size_t read_length) {
  if (!out_ || !*out_) return;

  /* Unmapped: * for the target columns; keep the query span and the best chain's stats */
  if (!result.mapped()) {
    std::stringstream ss;
    if (!result.mappings.empty()) {
      const auto& best = result.mappings[0];
      const auto s = tags::chain_span(best);
      ss << read_id << '\t' << read_length << '\t' << s.q_lo << '\t' << s.q_hi
         << "\t*\t*\t*\t*\t*\t0\t0\t" << best.mapq << "\tpn:Z:*";
      tags::write_chain(ss, 'U', best, s);
    } else {
      ss << read_id << '\t' << read_length << "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\tpn:Z:*\ttp:A:U";
    }
    tags::write_read(ss, result);
    *out_ << ss.str() << '\n';
    return;
  }

  if (result.mappings.empty()) return;

  std::size_t max_out = result.mappings.size();
  if (config_.primary_only) {
    max_out = 1;
  } else {
    max_out = std::min(max_out, config_.max_secondary + 1);
  }

  for (std::size_t i = 0; i < max_out; ++i) {
    const auto& mapping = result.mappings[i];
    const bool is_primary = (i == 0);
    const auto s = tags::chain_span(mapping);
    const bool on_path = mapping.coord_space != mapping::CoordSpace::kCanonical &&
                         mapping.path_id < graph_.pathCount();

    std::string walk;
    std::string path_name = "*";
    std::int64_t walk_len = 0, start = s.r_lo, end = s.r_hi;
    char strand = '+';
    constexpr std::string_view kRevSuffix = "_reverse";
    auto strip_reverse = [&](std::string& name) {
      const bool rev = name.size() > kRevSuffix.size() &&
                       name.compare(name.size() - kRevSuffix.size(), kRevSuffix.size(), kRevSuffix) == 0;
      if (rev) name.resize(name.size() - kRevSuffix.size());
      return rev;
    };
    const auto pid = static_cast<std::uint32_t>(mapping.path_id);
    if (on_path && !config_.anchor_nodes && graph_.pathStepCount(pid) > 0) {
      // Full walk: every path step from the step that holds the chain start to the step
      // that holds its end. Columns 7 to 9 are the walk length and the chain start and
      // end on the walk, as the GAF spec defines them. A "<name>_reverse" path gives a
      // reverse walk on strand '+'.
      const auto& starts = stepStarts(pid);
      const auto path_len = static_cast<std::int64_t>(starts.back());
      const std::int64_t lo = std::max(std::int64_t{0}, s.r_lo);
      const std::int64_t hi = std::max(lo + 1, std::min(s.r_hi, path_len));
      const auto last = starts.end() - 1;
      const std::size_t first_step = std::upper_bound(starts.begin(), last, lo) - starts.begin() - 1;
      const std::size_t last_step = std::upper_bound(starts.begin(), last, hi - 1) - starts.begin() - 1;
      const auto* steps = graph_.pathStepsBegin(pid);
      for (std::size_t k = first_step; k <= last_step; ++k) {
        const auto node = steps[k];
        walk += graph_.isReverse(node) ? '<' : '>';
        const auto name = graph_.name(node);
        walk.append(name.data(), name.size());
      }
      walk_len = static_cast<std::int64_t>(starts[last_step + 1] - starts[first_step]);
      start = lo - static_cast<std::int64_t>(starts[first_step]);
      end = hi - static_cast<std::int64_t>(starts[first_step]);
      path_name = std::string(graph_.pathName(pid));
      strip_reverse(path_name);
    } else if (on_path) {
      // Anchored nodes only (--gaf-anchor-nodes, or a graph without path steps): column 6
      // lists the nodes the anchors fall on, and columns 7 to 9 are the path length and
      // the positions on the path (forward strand), as in PAF.
      walk = buildPathString(mapping.anchors);
      path_name = std::string(graph_.pathName(pid));
      walk_len = static_cast<std::int64_t>(graph_.pathLength(pid));
      start = std::max(std::int64_t{0}, s.r_lo);
      end = walk_len > 0 ? std::min(s.r_hi, walk_len) : s.r_hi;
      if (strip_reverse(path_name)) {
        strand = '-';
        const auto flipped_start = walk_len - end;
        end = walk_len - start;
        start = flipped_start;
      }
    } else {
      // 1D canonical coordinates (sort chainer): no path. Column 6 lists the anchored
      // nodes; columns 8 and 9 are 1D positions.
      walk = buildPathString(mapping.anchors);
    }
    const std::int64_t block = end - start;

    std::stringstream ss;
    ss << read_id << '\t' << read_length << '\t' << s.q_lo << '\t' << s.q_hi << '\t' << strand << '\t' << walk
       << '\t' << walk_len << '\t' << start << '\t' << end << '\t' << tags::matches(mapping, block)
       << '\t' << block << '\t' << mapping.mapq << "\tpn:Z:" << path_name;
    tags::write_chain(ss, is_primary ? 'P' : 'S', mapping, s);

    // 1D layout interval in forward-node space (ls:f start, le:f end, lc:i component).
    // Reverse nodes map back to their forward partner with the offset flipped.
    if (config_.node_1d_coords && config_.component_ids) {
      const auto& coords = *config_.node_1d_coords;
      const auto& comp = *config_.component_ids;
      float ls = std::numeric_limits<float>::max();
      float le = std::numeric_limits<float>::lowest();
      for (const auto& a : mapping.anchors) {
        std::uint32_t fwd_id = a.node_id & ~1u;
        float fwd_start = coords[fwd_id];
        float fwd_end = coords[fwd_id + 1];
        float canon_span = fwd_end - fwd_start;
        float base_len = static_cast<float>(graph_.seqLen(fwd_id));
        float pos, pos_end;
        if (base_len < 1.0f) {
          pos = pos_end = fwd_start;
        } else if (a.node_id & 1u) {
          float frac = static_cast<float>(a.offset) / base_len;
          float frac_end = static_cast<float>(a.offset + a.length) / base_len;
          pos = fwd_start + (1.0f - frac_end) * canon_span;
          pos_end = fwd_start + (1.0f - frac) * canon_span;
        } else {
          float frac = static_cast<float>(a.offset) / base_len;
          float frac_end = static_cast<float>(a.offset + a.length) / base_len;
          pos = fwd_start + frac * canon_span;
          pos_end = fwd_start + frac_end * canon_span;
        }
        ls = std::min(ls, std::min(pos, pos_end));
        le = std::max(le, std::max(pos, pos_end));
      }
      tags::put_f(ss, "ls:f:%.1f", ls);
      tags::put_f(ss, "le:f:%.1f", le);
      ss << "\tlc:i:" << comp[mapping.anchors[0].node_id & ~1u];
    }

    if (is_primary) tags::write_read(ss, result);
    *out_ << ss.str() << '\n';
  }
}

const std::vector<std::uint64_t>& GafWriter::stepStarts(std::uint32_t path_id) const {
  auto it = step_starts_.find(path_id);
  if (it != step_starts_.end()) return it->second;
  const std::size_t n = graph_.pathStepCount(path_id);
  const auto* steps = graph_.pathStepsBegin(path_id);
  std::vector<std::uint64_t> starts(n + 1, 0);
  for (std::size_t k = 0; k < n; ++k) starts[k + 1] = starts[k] + graph_.seqLen(steps[k]);
  return step_starts_.emplace(path_id, std::move(starts)).first->second;
}

std::string GafWriter::buildPathString(const std::vector<mapping::ChainedAnchor>& anchors) const {
  if (anchors.empty()) return "*";

  std::string path;
  std::size_t last_node_id = std::numeric_limits<std::size_t>::max();

  for (const auto& anchor : anchors) {
    if (anchor.node_id == last_node_id) continue;

    auto node_name = graph_.name(anchor.node_id);
    char direction = graph_.isReverse(anchor.node_id) ? '<' : '>';

    path += direction;
    path.append(node_name.data(), node_name.size());
    last_node_id = anchor.node_id;
  }

  return path.empty() ? "*" : path;
}

}  // namespace panomap::io
