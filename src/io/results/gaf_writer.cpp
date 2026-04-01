/**
 * gaf_writer.cpp
 *
 * Writes mapping results directly to GAF format. Handles coordinate
 * transforms, path lookups, and tag generation internally.
 *
 * SPDX-License-Identifier: MIT
 */

#include "io/results/gaf_writer.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string_view>
#include <unordered_map>

#include "util/logging.hpp"

namespace piru::io {

GafWriter::GafWriter(const std::string& path, const index::FlatGraph& graph,
                     GafWriterConfig config)
    : out_(path, std::ios::out | std::ios::trunc), graph_(graph), config_(std::move(config)) {
  if (!out_) {
    LOG_ERROR("Failed to open GAF output: " + path);
  } else {
    out_ << "#query_name\tquery_len\tquery_start\tquery_end\tstrand\t"
         << "graph_path\tpath_len\tpath_start\tpath_end\t"
         << "matches\tblock_len\tmapq\t[optional_tags: pn:Z:path_name ...]\n";
  }
}

GafWriter::~GafWriter() { out_.flush(); }

void GafWriter::write(const mapping::ReadMapResult& result,
                      const std::string& read_id,
                      std::size_t read_length) {
  if (!out_ || result.mappings.empty()) return;

  /* Determine how many mappings to output */
  std::size_t max_out = result.mappings.size();
  if (config_.primary_only) {
    max_out = 1;
  } else {
    max_out = std::min(max_out, config_.max_secondary + 1);
  }

  const double primary_score = result.mappings[0].chain_score;
  const double min_score = primary_score * config_.min_secondary_ratio;

  for (std::size_t i = 0; i < max_out; ++i) {
    const auto& mapping = result.mappings[i];
    bool is_primary = (i == 0);

    // Skip weak secondaries
    if (!is_primary && mapping.chain_score < min_score) continue;

    bool is_canonical = (mapping.coord_space == mapping::CoordSpace::kCanonical);

    /* Query span */
    std::size_t min_query = mapping.anchors.front().read_pos;
    std::size_t max_query = mapping.anchors.front().read_pos;
    for (const auto& a : mapping.anchors) {
      min_query = std::min(min_query, static_cast<std::size_t>(a.read_pos));
      max_query = std::max(max_query, static_cast<std::size_t>(a.read_pos) + a.length);
    }

    /* Target info + ref span */
    std::string target_path;
    std::size_t target_length = 0;
    std::int64_t min_ref, max_ref;
    char strand = '+';

    if (is_canonical) {
      target_path = "*";
      min_ref = mapping.anchors.front().ref_coord;
      max_ref = mapping.anchors.front().ref_coord;
      for (const auto& a : mapping.anchors) {
        min_ref = std::min(min_ref, a.ref_coord);
        max_ref = std::max(max_ref, a.ref_coord + static_cast<std::int64_t>(a.length));
      }
      // TODO: pass canonical 1D span through config
      target_length = 0;
    } else {
      target_path = getPathName(mapping.path_id);
      target_length = getPathLength(mapping.path_id);
      auto path_len = static_cast<std::int64_t>(target_length);

      min_ref = std::max(std::int64_t{0}, mapping.anchors.front().ref_coord);
      max_ref = mapping.anchors.front().ref_coord;
      for (const auto& a : mapping.anchors) {
        min_ref = std::min(min_ref, std::max(std::int64_t{0}, a.ref_coord));
        auto end = a.ref_coord + static_cast<std::int64_t>(a.length);
        max_ref = std::max(max_ref, std::min(end, path_len));
      }

      // Strand flipping for reverse paths
      constexpr std::string_view kRevSuffix = "_reverse";
      if (target_path.size() > kRevSuffix.size() &&
          target_path.substr(target_path.size() - kRevSuffix.size()) == kRevSuffix) {
        strand = '-';
        target_path = target_path.substr(0, target_path.size() - kRevSuffix.size());
        auto flipped_start = path_len - max_ref;
        auto flipped_end = path_len - min_ref;
        min_ref = flipped_start;
        max_ref = flipped_end;
      }
    }

    /* Matches / block length from anchor coverage */
    std::uint64_t total_anchor_len = 0;
    for (const auto& a : mapping.anchors) total_anchor_len += a.length;

    /* Graph path string (node walk) */
    std::string graph_path = buildPathString(mapping.anchors);

    /* Build GAF line */
    std::stringstream ss;
    const std::string& path_col = graph_path.empty() ? target_path : graph_path;
    ss << read_id << '\t' << read_length << '\t' << min_query << '\t' << max_query
       << '\t' << strand << '\t' << path_col << '\t' << target_length
       << '\t' << min_ref << '\t' << max_ref
       << '\t' << total_anchor_len << '\t' << (max_ref - min_ref) << '\t' << mapping.mapq;

    // pn:Z: path name tag (when graph path is used in col 6)
    if (!graph_path.empty() && !target_path.empty() && target_path != "*") {
      ss << "\tpn:Z:" << target_path;
    }

    // tp:A: tag
    ss << "\ttp:A:" << (is_primary ? 'P' : 'S');

    // cs:i: chain score
    ss << "\tcs:i:" << static_cast<int>(mapping.chain_score);

    // an:i: anchor count
    ss << "\tan:i:" << mapping.anchors.size();

    // Per-read timing (primary only)
    if (is_primary) {
      ss << "\tck:i:" << result.chunks_processed;
      char dt_buf[32];
      std::snprintf(dt_buf, sizeof(dt_buf), "dt:f:%.6f", result.processing_time_sec);
      ss << '\t' << dt_buf;
    }

    // ROI tags (primary only)
    if (is_primary && result.roi_overlap >= 0.0) {
      ss << "\tro:f:" << result.roi_overlap;
      ss << "\trd:A:" << (result.roi_keep ? 'K' : 'R');
    }


    out_ << ss.str() << '\n';
  }
}

std::string GafWriter::buildPathString(
    const std::vector<mapping::ChainedAnchor>& anchors) const {
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

std::string GafWriter::getPathName(std::size_t path_id) const {
  if (path_id >= graph_.pathCount()) return "unknown";
  return std::string(graph_.pathName(path_id));
}

std::size_t GafWriter::getPathLength(std::size_t path_id) const {
  if (path_id < graph_.pathCount() && graph_.pathLength(static_cast<std::uint32_t>(path_id)) > 0) {
    return graph_.pathLength(static_cast<std::uint32_t>(path_id));
  }

  if (!path_lengths_computed_) {
    path_lengths_.resize(graph_.pathCount());

    for (std::uint32_t p = 0; p < graph_.pathCount(); ++p) {
      if (graph_.pathLength(p) > 0) {
        path_lengths_[p] = graph_.pathLength(p);
        continue;
      }
      // Compute from step node lengths (no overlap tracking in FlatGraph)
      std::size_t len = 0;
      const auto* steps = graph_.pathStepsBegin(p);
      std::size_t num_steps = graph_.pathStepCount(p);
      for (std::size_t s = 0; s < num_steps; ++s) {
        len += graph_.seqLen(steps[s]);
      }
      path_lengths_[p] = len;
    }
    path_lengths_computed_ = true;
  }

  return (path_id < path_lengths_.size()) ? path_lengths_[path_id] : 0;
}

}  // namespace piru::io
