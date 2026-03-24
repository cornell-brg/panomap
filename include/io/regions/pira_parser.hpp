// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

namespace piru::io {

/// A single annotated region with optional 1D interval.
struct PiraRegion {
  std::string label;                       // e.g. "path:start-end"
  double coord_1d_start{0.0};             // 1D canonical interval start (* = not set)
  double coord_1d_end{0.0};               // 1D canonical interval end
  bool has_1d{false};                      // whether 1D coords are present
  std::unordered_set<std::size_t> nodes;   // ROI node IDs
};

/// Parsed .pira file contents.
struct PiraFile {
  std::vector<PiraRegion> regions;
  std::unordered_set<std::size_t> all_nodes;  // union of all region nodes
};

/// Parse a .pira annotation file (v1 or v2 format).
/// v1: "label node1,node2,..."
/// v2: "label\t1d_start\t1d_end\tnode1,node2,..."
PiraFile parse_pira_v2(const std::string& filepath);

/// Legacy: parse and return flat node set (v1 compat).
std::unordered_set<std::size_t> parse_pira(const std::string& filepath);

}  // namespace piru::io
