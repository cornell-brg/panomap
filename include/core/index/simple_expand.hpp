// SPDX-License-Identifier: MIT
// Simple +-expand transformation: ImportedGraph -> FlatGraph with exactly 2x nodes.

#pragma once

#include "core/index/flat_graph.hpp"
#include "core/io/graphs/graph.hpp"

namespace piru::index {

// Expand ImportedGraph to FlatGraph.
// Each original node becomes 2 nodes (forward + reverse).
// Node ID scheme: original_index * 2 = forward, original_index * 2 + 1 = reverse.
// Paths are duplicated: forward + reverse for each original path.
FlatGraph simpleExpandFlat(const piru::io::ImportedGraph& imported);

// Helper functions for node ID arithmetic
inline std::size_t forwardNodeId(std::size_t original_index) { return original_index * 2; }
inline std::size_t reverseNodeId(std::size_t original_index) { return original_index * 2 + 1; }
inline std::size_t originalIndex(std::size_t aln_id) { return aln_id / 2; }
inline bool isReverseNode(std::size_t aln_id) { return (aln_id & 1) == 1; }

}  // namespace piru::index
