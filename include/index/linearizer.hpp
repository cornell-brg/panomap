// SPDX-License-Identifier: MIT
// Linear coordinate type for graph node linearization.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace piru::index {

// Linear coordinate assigned to a node.
// A node may have multiple coordinates if it appears on multiple paths
// or appears multiple times on the same path (cycles).
struct LinearCoordinate {
  std::size_t path_id;     // Which reference path
  std::int64_t ref_coord;  // Linear position along that path

  LinearCoordinate(std::size_t pid, std::int64_t coord) : path_id(pid), ref_coord(coord) {}
};

// Transposed coordinate index: [path_id][node_id] -> list of ref_coords.
// A node can appear multiple times on the same path (cycles/repeats).
// Built from the per-node coords for efficient per-path lookups.
using TransposedCoords = std::vector<std::vector<std::vector<std::int64_t>>>;

// Build transposed view from per-node coords.
// Input:  coords[node_id] -> vector<LinearCoordinate{path_id, ref_coord}>
// Output: result[path_id][node_id] -> vector<ref_coord>
inline TransposedCoords buildTransposedCoords(
    const std::vector<std::vector<LinearCoordinate>>& coords) {
  // Find max path_id
  std::size_t max_path = 0;
  for (const auto& node_coords : coords) {
    for (const auto& lc : node_coords) {
      if (lc.path_id > max_path) max_path = lc.path_id;
    }
  }

  TransposedCoords result(max_path + 1, std::vector<std::vector<std::int64_t>>(coords.size()));

  for (std::size_t node_id = 0; node_id < coords.size(); ++node_id) {
    for (const auto& lc : coords[node_id]) {
      result[lc.path_id][node_id].push_back(lc.ref_coord);
    }
  }

  return result;
}

}  // namespace piru::index
