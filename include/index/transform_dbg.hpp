// SPDX-License-Identifier: MIT
// Transform ImportedGraph (DBG) into directional AlnGraph.

#pragma once

#include <cstddef>

#include "index/aln_graph.hpp"
#include "io/graphs/graph.hpp"

namespace piru::index {

// Transform a DBG-style ImportedGraph into an AlnGraph with directional nodes.
// Trims node sequences by k_delta = (graph_k - pore_k) to preserve only pore_k-1 overlap.
// Adjusts edge overlaps accordingly. Does not perform pseudo-linearization.
//
// Params:
//   imported: Input DBG graph (assumes k-1 overlaps between nodes)
//   graph_k: DBG k-mer size (e.g., 63 for Bifrost)
//   pore_k: Pore model k-mer size (e.g., 5 for R9.4)
//
// Throws: std::invalid_argument if pore_k > graph_k
AlnGraph transformDbg(const io::ImportedGraph& imported, std::size_t graph_k, std::size_t pore_k);

}  // namespace piru::index
