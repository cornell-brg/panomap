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
#include <iostream>
#include <limits>
#include <sstream>
#include <string_view>
#include <unordered_map>

#include "util/logging.hpp"

namespace piru::io {

namespace {
void write_header(std::ostream& out) {
  out << "#query_name\tquery_len\tquery_start\tquery_end\tstrand\t"
      << "graph_path\tpath_len\tpath_start\tpath_end\t"
      << "matches\tblock_len\tmapq\t[optional_tags: pn:Z:path_name ...]\n";
}
}  // namespace

GafWriter::GafWriter(const std::string& path, const index::FlatGraph& graph, GafWriterConfig config)
    : file_(path, std::ios::out | std::ios::trunc), out_(&file_), graph_(graph), config_(std::move(config)) {
  if (!file_) {
    LOG_ERROR("Failed to open GAF output: " + path);
  } else {
    write_header(*out_);
  }
}

GafWriter::GafWriter(const index::FlatGraph& graph, GafWriterConfig config)
    : out_(&std::cout), graph_(graph), config_(std::move(config)) {
  write_header(*out_);
}

GafWriter::~GafWriter() { out_->flush(); }

void GafWriter::write(const mapping::ReadMapResult& result, const std::string& read_id,
                      std::size_t read_length) {
  if (!out_ || !*out_) return;

  /* Unmapped: output with * ref coordinates but keep query span and chain stats */
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
         << "\t*\t*\t*\t*\t*\t0\t0\t" << best.mapq;
      ss << "\tpn:Z:*\ttp:A:U";
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
      ss << read_id << '\t' << read_length << "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0";
      ss << "\tpn:Z:*\ttp:A:U";
    }
    // ws:f: weighted standout score from mapping decision
    // nc:i: number of competitive chains (after secondary ratio filter)
    // df:A: decision path: G=gate, S=standout, F=fallback, U=unmapped
    {
      char ws_buf[32];
      std::snprintf(ws_buf, sizeof(ws_buf), "ws:f:%.3f", result.standout);
      ss << '\t' << ws_buf;
      ss << "\tnc:i:" << result.mappings.size();
      ss << "\tdf:A:" << static_cast<char>(result.decision_path);
    }
    ss << "\tck:i:" << result.chunks_processed;
    {
      char dt_buf[32];
      std::snprintf(dt_buf, sizeof(dt_buf), "dt:f:%.6f", result.processing_time_sec);
      ss << '\t' << dt_buf;
    }
    *out_ << ss.str() << '\n';
    return;
  }

  if (result.mappings.empty()) return;

  /* Determine how many mappings to output */
  std::size_t max_out = result.mappings.size();
  if (config_.primary_only) {
    max_out = 1;
  } else {
    max_out = std::min(max_out, config_.max_secondary + 1);
  }

  for (std::size_t i = 0; i < max_out; ++i) {
    const auto& mapping = result.mappings[i];
    bool is_primary = (i == 0);

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
    ss << read_id << '\t' << read_length << '\t' << min_query << '\t' << max_query << '\t' << strand
       << '\t' << path_col << '\t' << target_length << '\t' << min_ref << '\t' << max_ref << '\t'
       << total_anchor_len << '\t' << (max_ref - min_ref) << '\t' << mapping.mapq;

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

    // se:f: score per event span (score / query_span)
    {
      auto qspan = max_query - min_query;
      if (qspan > 0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "se:f:%.3f", mapping.chain_score / static_cast<double>(qspan));
        ss << '\t' << buf;
      }
    }

    // ad:f: anchor density (anchors / ref_span)
    {
      auto rspan = max_ref - min_ref;
      if (rspan > 0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "ad:f:%.4f",
                      static_cast<double>(mapping.anchors.size()) / static_cast<double>(rspan));
        ss << '\t' << buf;
      }
    }

    // Canonical 1D interval in forward-node space (ci:f: start, ce:f: end, cc:i: component_id).
    // Reverse nodes map back to their forward partner with offset flipping.
    if (config_.node_1d_coords && config_.component_ids) {
      const auto& coords = *config_.node_1d_coords;
      const auto& comp = *config_.component_ids;
      float ci = std::numeric_limits<float>::max();
      float ce = std::numeric_limits<float>::lowest();
      for (const auto& a : mapping.anchors) {
        std::uint32_t fwd_id = a.node_id & ~1u;
        float fwd_start = coords[fwd_id];
        float fwd_end = coords[fwd_id + 1];
        float canon_span = fwd_end - fwd_start;
        float base_len = static_cast<float>(graph_.seqLen(fwd_id));
        // Interpolate offset within canonical span using base length ratio
        float pos, pos_end;
        if (base_len < 1.0f) {
          pos = pos_end = fwd_start;
        } else if (a.node_id & 1u) {
          // Reverse node: offset 0 on rev = fwd_end, increasing offset moves toward fwd_start
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
        ci = std::min(ci, std::min(pos, pos_end));
        ce = std::max(ce, std::max(pos, pos_end));
      }
      std::uint32_t cc = comp[mapping.anchors[0].node_id & ~1u];
      char ci_buf[32], ce_buf[32];
      std::snprintf(ci_buf, sizeof(ci_buf), "ci:f:%.1f", ci);
      std::snprintf(ce_buf, sizeof(ce_buf), "ce:f:%.1f", ce);
      ss << '\t' << ci_buf << '\t' << ce_buf << "\tcc:i:" << cc;
    }

    // Per-read timing (primary only)
    if (is_primary) {
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



    *out_ << ss.str() << '\n';
  }
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
