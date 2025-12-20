// SPDX-License-Identifier: MIT
// Linearization interface for assigning linear coordinates to graph nodes.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"

namespace piru::index {

// Linear coordinate assigned to a node.
// A node may have multiple coordinates if it appears on multiple paths
// or appears multiple times on the same path (cycles).
struct LinearCoordinate {
  std::size_t path_id;   // Which reference path (or chain ID for superbubble backend)
  std::int64_t ref_coord; // Linear position along that path/chain

  LinearCoordinate(std::size_t pid, std::int64_t coord) : path_id(pid), ref_coord(coord) {}
};

// Interface for assigning linear coordinates to graph nodes.
// Different backends implement different linearization strategies:
// - SuperbubbleLinearizer: Uses pseudo-linearization (chain IDs and local coordinates)
// - PathWalkLinearizer: Walks reference paths to assign global coordinates
class Linearizer {
public:
  virtual ~Linearizer() = default;

  // Assign linear coordinates to all nodes in the graph using signal sizes.
  // signal_sizes[node_id] = number of signal samples for that node.
  // Returns a vector of coordinates per node (indexed by node ID).
  // Nodes may have zero coordinates (not on any path/chain),
  // one coordinate (on single path/chain), or multiple coordinates
  // (on multiple paths or multiple occurrences on same path).
  virtual std::vector<std::vector<LinearCoordinate>> linearize(
      const AlnGraph& graph, const std::vector<std::size_t>& signal_sizes) const = 0;

  // Backend name for logging and debugging.
  virtual std::string name() const = 0;
};

}  // namespace piru::index
