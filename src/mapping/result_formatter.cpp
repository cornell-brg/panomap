// SPDX-License-Identifier: MIT

#include "mapping/result_formatter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace piru::mapping {

ResultFormatter::ResultFormatter(const index::AlnGraph& graph,
                                 ResultFormatterConfig config)
    : graph_(graph), config_(std::move(config)) {}

std::vector<io::AlignmentResult> ResultFormatter::format(
    const ReadMapResult& map_result,
    const std::string& read_id,
    std::size_t read_length) const {

  std::vector<io::AlignmentResult> results;

  if (map_result.mappings.empty()) {
    return results;
  }

  // Determine how many mappings to output
  std::size_t num_mappings = map_result.mappings.size();
  if (config_.primary_only) {
    num_mappings = 1;
  } else {
    num_mappings = std::min(num_mappings, config_.max_secondary + 1);
  }

  results.reserve(num_mappings);

  // Get primary and secondary scores for MAPQ calculation
  double primary_score = map_result.mappings[0].chain_score;
  double secondary_score = (map_result.mappings.size() > 1)
                               ? map_result.mappings[1].chain_score
                               : 0.0;

  for (std::size_t i = 0; i < num_mappings; ++i) {
    bool is_primary = (i == 0);
    const auto& mapping = map_result.mappings[i];

    // Format mapping to result
    results.push_back(formatMapping(mapping, read_id, read_length, is_primary));

    // Set MAPQ based on primary vs secondary score gap
    if (is_primary) {
      results.back().mapq = estimateMapQ(primary_score, secondary_score);
    } else {
      // Secondary alignments get lower MAPQ
      results.back().mapq = 0;
    }

    // Add tp:A:P (primary) or tp:A:S (secondary) tag
    results.back().optional_fields.push_back(
        is_primary ? "tp:A:P" : "tp:A:S");

    // Add chain score tag
    results.back().optional_fields.push_back(
        "cs:i:" + std::to_string(static_cast<int>(mapping.chain_score)));

    // Add anchor count tag
    results.back().optional_fields.push_back(
        "an:i:" + std::to_string(mapping.anchors.size()));

    // Add alignment score tag if available
    if (results.back().alignment_score.has_value()) {
      results.back().optional_fields.push_back(
          "as:f:" + std::to_string(results.back().alignment_score.value()));
    }

    // Set graph path traversal (for GAF column 6)
    results.back().graph_path = buildPathString(mapping.anchors);
  }

  return results;
}

io::AlignmentResult ResultFormatter::formatMapping(
    const Mapping& mapping,
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
    min_query = std::min(min_query, anchor.read_pos);
    max_query = std::max(max_query, anchor.read_pos + anchor.target.length);
  }
  result.query_start = min_query;
  result.query_end = max_query;

  // Target info - use path from first anchor (all anchors in chain should be same path)
  std::size_t path_id = mapping.anchors.front().path_id;
  result.target_path = getPathName(path_id);
  result.target_length = computePathLength(path_id);

  // Find reference span from anchors
  std::int64_t min_ref = mapping.anchors.front().ref_coord;
  std::int64_t max_ref = mapping.anchors.front().ref_coord;
  for (const auto& anchor : mapping.anchors) {
    min_ref = std::min(min_ref, anchor.ref_coord);
    max_ref = std::max(max_ref, anchor.ref_coord + static_cast<std::int64_t>(anchor.target.length));
  }
  result.target_start = static_cast<std::uint64_t>(min_ref);
  result.target_end = static_cast<std::uint64_t>(max_ref);

  // Strand - determine from path direction
  // For now, assume forward (+). Path names ending with "_reverse" indicate reverse.
  result.strand = '+';
  if (result.target_path.size() > 8 &&
      result.target_path.substr(result.target_path.size() - 8) == "_reverse") {
    result.strand = '-';
  }

  // Approximate matches and block length from anchor coverage
  // matches = sum of anchor lengths (approximate)
  std::uint64_t total_anchor_length = 0;
  for (const auto& anchor : mapping.anchors) {
    total_anchor_length += anchor.target.length;
  }
  result.matches = total_anchor_length;
  result.alignment_block_length = result.target_end - result.target_start;

  // Build mappings from anchors (for GAF path output)
  result.mappings.reserve(mapping.anchors.size());
  for (const auto& anchor : mapping.anchors) {
    io::AlignmentResult::Mapping m;
    m.node_id = anchor.target.node_id;
    m.offset = static_cast<std::uint32_t>(anchor.target.offset);
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

std::string ResultFormatter::buildPathString(
    const std::vector<SeedAnchor>& anchors) const {

  if (anchors.empty()) {
    return "*";
  }

  std::string path;
  std::size_t last_node_id = std::numeric_limits<std::size_t>::max();

  for (const auto& anchor : anchors) {
    // Skip consecutive anchors on the same node
    if (anchor.target.node_id == last_node_id) {
      continue;
    }

    // Get original node ID from AlnGraph
    const auto& node = graph_.node(anchor.target.node_id);
    const std::string& original_id = node.original_id.empty()
                                         ? node.label
                                         : node.original_id;

    // Determine direction (> for forward, < for reverse)
    char direction = node.is_reverse ? '<' : '>';

    path += direction;
    path += original_id;

    last_node_id = anchor.target.node_id;
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
      std::size_t cumulative_len = 0;

      for (std::size_t step_idx = 0; step_idx < path.steps.size(); ++step_idx) {
        const auto& step = path.steps[step_idx];

        auto it = label_to_idx.find(step.node_id);
        if (it == label_to_idx.end()) {
          continue;
        }

        const auto& node = graph_.node(it->second);
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

int ResultFormatter::estimateMapQ(double primary_score, double secondary_score) const {
  // Simple MAPQ estimation based on score gap
  // Similar to minimap2's approach: mapq = 40 * (1 - secondary/primary)
  // Clamped to [0, 60]

  if (primary_score <= 0) {
    return 0;
  }

  if (secondary_score <= 0) {
    // No secondary alignment - high confidence
    return 60;
  }

  double ratio = secondary_score / primary_score;
  int mapq = static_cast<int>(40.0 * (1.0 - ratio));

  return std::clamp(mapq, 0, 60);
}

}  // namespace piru::mapping
