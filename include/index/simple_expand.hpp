// SPDX-License-Identifier: MIT
// Simple ±expand transformation: ImportedGraph → AlnGraph with exactly 2x nodes.

#pragma once

#include "index/aln_graph.hpp"
#include "io/graphs/graph.hpp"

namespace piru::index {

// Simple expansion of ImportedGraph to AlnGraph.
// Each original node becomes exactly 2 AlnGraph nodes (forward + reverse).
// Node ID scheme: original_index * 2 = forward, original_index * 2 + 1 = reverse.
// Paths are duplicated: forward path + reverse path for each original path.
AlnGraph simpleExpand(const piru::io::ImportedGraph& imported);

// Helper functions for node ID arithmetic
inline std::size_t forwardNodeId(std::size_t original_index) { return original_index * 2; }
inline std::size_t reverseNodeId(std::size_t original_index) { return original_index * 2 + 1; }
inline std::size_t originalIndex(std::size_t aln_id) { return aln_id / 2; }
inline bool isReverseNode(std::size_t aln_id) { return (aln_id & 1) == 1; }

}  // namespace piru::index
