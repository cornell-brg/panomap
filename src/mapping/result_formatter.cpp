// SPDX-License-Identifier: MIT

#include "mapping/result_formatter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>
#include <unordered_map>

#include "util/logging.hpp"

namespace piru::mapping {

ResultFormatter::ResultFormatter(const index::AlnGraph& graph, ResultFormatterConfig config)
    : graph_(graph), config_(std::move(config)) {}

std::vector<io::AlignmentResult> ResultFormatter::format(const ReadMapResult& map_result,
                                                         const std::string& read_id,
                                                         std::size_t read_length) const {
  std::vector<io::AlignmentResult> results;

  if (map_result.mappings.empty()) {
    return results;
  }

  // Determine how many mappings to output
  std::size_t max_mappings = map_result.mappings.size();
  if (config_.primary_only) {
    max_mappings = 1;
  } else {
    max_mappings = std::min(max_mappings, config_.max_secondary + 1);
  }

  // Primary chain score for filtering secondaries
  const double primary_score = map_result.mappings[0].chain_score;
  const double min_score = primary_score * config_.min_secondary_ratio;

  results.reserve(max_mappings);

  for (std::size_t i = 0; i < max_mappings; ++i) {
    const auto& mapping = map_result.mappings[i];
    bool is_primary = (i == 0);

    // Skip secondaries below the score threshold
    if (!is_primary && mapping.chain_score < min_score) {
      continue;
    }

    // Format mapping to result
    results.push_back(formatMapping(mapping, read_id, read_length, is_primary));

    // Set MAPQ based on chain score
    results.back().mapq = estimateMapQ(mapping.chain_score, 0);

    // Add tp:A:P (primary) or tp:A:S (secondary) tag
    results.back().optional_fields.push_back(is_primary ? "tp:A:P" : "tp:A:S");

    // Add chain score tag
    results.back().optional_fields.push_back("cs:i:" +
                                             std::to_string(static_cast<int>(mapping.chain_score)));

    // Add anchor count tag
    results.back().optional_fields.push_back("an:i:" + std::to_string(mapping.anchors.size()));

    // Add alignment score tag if available
    if (results.back().alignment_score.has_value()) {
      results.back().optional_fields.push_back(
          "as:f:" + std::to_string(results.back().alignment_score.value()));
    }

    // Set graph path traversal (for GAF column 6)
    results.back().graph_path = buildPathString(mapping.anchors);

    // ROI classification tags (primary only)
    if (is_primary && map_result.roi_overlap >= 0.0) {
      results.back().optional_fields.push_back("ro:f:" + std::to_string(map_result.roi_overlap));
      results.back().optional_fields.push_back(std::string("rd:A:") +
                                               (map_result.roi_keep ? "K" : "R"));
    }
  }

  return results;
}

io::AlignmentResult ResultFormatter::formatMapping(const Mapping& mapping,
                                                   const std::string& read_id,
                                                   std::size_t read_length,
                                                   bool /*is_primary*/) const {
  io::AlignmentResult result;

  if (mapping.anchors.empty()) {
    return result;
  }

  // Query info
  result.query_name = read_id;
  result.query_length = read_length;

  // Find query span from anchors
  std::size_t min_query = mapping.anchors.front().read_pos;
  std::size_t max_query = mapping.anchors.front().read_pos;
  for (const auto& anchor : mapping.anchors) {
    min_query = std::min(min_query, static_cast<std::size_t>(anchor.read_pos));
    max_query = std::max(max_query, static_cast<std::size_t>(anchor.read_pos) + anchor.length);
  }
  result.query_start = min_query;
  result.query_end = max_query;

  // Target info - use path from first anchor (all anchors in chain should be same path)
  std::size_t path_id = mapping.anchors.front().path_id;
  result.target_path = getPathName(path_id);
  result.target_length = computePathLength(path_id);

  // Find reference span from anchors
  // Clamp to path bounds (seeds near node boundaries can span past path end)
  const std::int64_t path_len = static_cast<std::int64_t>(result.target_length);
  std::int64_t min_ref = std::max(std::int64_t{0}, mapping.anchors.front().ref_coord);
  std::int64_t max_ref = mapping.anchors.front().ref_coord;
  for (const auto& anchor : mapping.anchors) {
    min_ref = std::min(min_ref, std::max(std::int64_t{0}, anchor.ref_coord));
    std::int64_t anchor_end = anchor.ref_coord + static_cast<std::int64_t>(anchor.length);
    max_ref = std::max(max_ref, std::min(anchor_end, path_len));
  }

  // Strand and coordinate handling
  // Path names ending with "_reverse" indicate reverse complement orientation.
  // For reverse paths: strip suffix and flip coordinates to original path space.
  result.strand = '+';
  constexpr std::string_view kReverseSuffix = "_reverse";
  if (result.target_path.size() > kReverseSuffix.size() &&
      result.target_path.substr(result.target_path.size() - kReverseSuffix.size()) ==
          kReverseSuffix) {
    result.strand = '-';

    // Strip "_reverse" suffix to get original path name
    result.target_path =
        result.target_path.substr(0, result.target_path.size() - kReverseSuffix.size());

    // Flip coordinates: convert from reverse path space to original path space
    // new_start = path_len - old_end
    // new_end = path_len - old_start
    std::int64_t path_len = static_cast<std::int64_t>(result.target_length);
    std::int64_t flipped_start = path_len - max_ref;
    std::int64_t flipped_end = path_len - min_ref;

    min_ref = flipped_start;
    max_ref = flipped_end;
  }

  result.target_start = static_cast<std::uint64_t>(min_ref);
  result.target_end = static_cast<std::uint64_t>(max_ref);

  // Approximate matches and block length from anchor coverage
  // matches = sum of anchor lengths (approximate)
  std::uint64_t total_anchor_length = 0;
  for (const auto& anchor : mapping.anchors) {
    total_anchor_length += anchor.length;
  }
  result.matches = total_anchor_length;
  result.alignment_block_length = result.target_end - result.target_start;

  // Build mappings from anchors (for GAF path output)
  result.mappings.reserve(mapping.anchors.size());
  for (const auto& anchor : mapping.anchors) {
    io::AlignmentResult::Mapping m;
    m.node_id = anchor.node_id;
    m.offset = static_cast<std::uint32_t>(anchor.offset);
    m.is_reverse = false;  // TODO: determine from path step direction
    // No edits for chain-only mode
    result.mappings.push_back(std::move(m));
  }

  // Add alignment score if available (already computed in process_read)
  if (mapping.alignment_cost.has_value()) {
    result.alignment_score = mapping.alignment_cost.value();
    if (mapping.normalized_alignment_cost.has_value()) {
      result.normalized_alignment_score = mapping.normalized_alignment_cost.value();
    }
  }

  return result;
}

std::string ResultFormatter::buildPathString(const std::vector<ChainedAnchor>& anchors) const {
  if (anchors.empty()) {
    return "*";
  }

  std::string path;
  std::size_t last_node_id = std::numeric_limits<std::size_t>::max();

  for (const auto& anchor : anchors) {
    // Skip consecutive anchors on the same node
    if (anchor.node_id == last_node_id) {
      continue;
    }

    // Get original node ID from AlnGraph
    const auto& node = graph_.node(anchor.node_id);
    const std::string& original_id = node.original_id.empty() ? node.label : node.original_id;

    // Determine direction (> for forward, < for reverse)
    char direction = node.is_reverse ? '<' : '>';

    path += direction;
    path += original_id;

    last_node_id = anchor.node_id;
  }

  return path.empty() ? "*" : path;
}

const std::string& ResultFormatter::getPathName(std::size_t path_id) const {
  static const std::string kUnknown = "unknown";

  if (path_id >= graph_.pathCount()) {
    return kUnknown;
  }

  return graph_.paths()[path_id].name;
}

std::size_t ResultFormatter::computePathLength(std::size_t path_id) const {
  // Use precomputed path length if available (simple pipeline)
  if (path_id < graph_.pathCount() && graph_.paths()[path_id].length > 0) {
    return graph_.paths()[path_id].length;
  }

  if (!path_lengths_computed_) {
    // Compute all path lengths on first access
    const auto& paths = graph_.paths();
    path_lengths_.resize(paths.size());

    // Build label -> node index mapping
    std::unordered_map<std::string, std::size_t> label_to_idx;
    for (std::size_t i = 0; i < graph_.nodeCount(); ++i) {
      label_to_idx[graph_.node(i).label] = i;
    }

    for (std::size_t p = 0; p < paths.size(); ++p) {
      const auto& path = paths[p];

      // Skip if already set
      if (path.length > 0) {
        path_lengths_[p] = path.length;
        continue;
      }

      std::size_t cumulative_len = 0;

      for (std::size_t step_idx = 0; step_idx < path.steps.size(); ++step_idx) {
        const auto& step = path.steps[step_idx];

        auto it = label_to_idx.find(step.node_id);
        if (it == label_to_idx.end()) {
          continue;
        }

        const std::size_t node_idx = it->second;

        // Sequence-based path length computation
        const auto& node = graph_.node(node_idx);
        std::size_t node_len = node.sequence.size();
        std::size_t overlap = 0;

        if (step_idx < path.steps.size() - 1 && step_idx < path.overlaps.size()) {
          overlap = path.overlaps[step_idx];
        }

        cumulative_len += node_len - overlap;
      }

      path_lengths_[p] = cumulative_len;
    }

    path_lengths_computed_ = true;
  }

  if (path_id >= path_lengths_.size()) {
    return 0;
  }

  return path_lengths_[path_id];
}

int ResultFormatter::estimateMapQ(double primary_score, double /*secondary_score*/) const {
  // Use raw chain score as MAPQ.
  // TODO: Refine MAPQ calculation when alignment scoring is added.
  return std::max(0, static_cast<int>(primary_score));
}

}  // namespace piru::mapping
