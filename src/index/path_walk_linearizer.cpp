// SPDX-License-Identifier: MIT
// Path-guided linearization implementation.

#include "index/path_walk_linearizer.hpp"

#include <stdexcept>
#include <unordered_map>

namespace piru::index {

std::vector<std::vector<LinearCoordinate>> PathWalkLinearizer::linearize(
    const AlnGraph& graph, const std::vector<std::size_t>& signal_sizes) const {
  const std::size_t n = graph.nodeCount();
  const auto& paths = graph.paths();

  // Error if no paths.
  if (paths.empty()) {
    throw std::runtime_error(
        "PathWalkLinearizer: graph has no reference paths. "
        "Use SuperbubbleLinearizer or provide a graph with paths.");
  }

  // Validate signal_sizes.
  if (signal_sizes.size() != n) {
    throw std::runtime_error("PathWalkLinearizer: signal_sizes.size() != nodeCount()");
  }

  // Build label → node index mapping.
  std::unordered_map<std::string, std::size_t> label_to_idx;
  for (std::size_t i = 0; i < n; ++i) {
    label_to_idx[graph.node(i).label] = i;
  }

  // Initialize result: each node can have multiple coordinates.
  std::vector<std::vector<LinearCoordinate>> result(n);

  // Walk each path and assign coordinates.
  for (std::size_t path_id = 0; path_id < paths.size(); ++path_id) {
    const auto& path = paths[path_id];
    std::int64_t cumulative_pos = 0;

    for (std::size_t step_idx = 0; step_idx < path.steps.size(); ++step_idx) {
      const auto& step = path.steps[step_idx];

      // Find node index from label.
      auto it = label_to_idx.find(step.node_id);
      if (it == label_to_idx.end()) {
        throw std::runtime_error("PathWalkLinearizer: path step references unknown node: " +
                                 step.node_id);
      }
      const std::size_t node_idx = it->second;

      // Record coordinate for this node occurrence.
      // Note: This is the start position of the node on this path.
      // When expanding seeds, we'll add the offset within the node.
      result[node_idx].emplace_back(path_id, cumulative_pos);

      // Advance cumulative position by signal size (no overlap math needed).
      cumulative_pos += static_cast<std::int64_t>(signal_sizes[node_idx]);
    }
  }

  return result;
}

}  // namespace piru::index
