// SPDX-License-Identifier: MIT
// Abstract interface for transforming ImportedGraph to AlnGraph, specifically for VG graphs.

#pragma once

#include "io/graphs/graph.hpp"
#include "index/aln_graph.hpp"

namespace piru::index {

// Structure to hold statistics about the graph transformation.
struct TransformStats {
    // Number of nodes in the original ImportedGraph.
    std::size_t original_node_count{0};
    // Number of nodes in the transformed AlnGraph.
    std::size_t transformed_node_count{0};
    // Number of edges in the original ImportedGraph.
    std::size_t original_edge_count{0};
    // Number of edges in the transformed AlnGraph.
    std::size_t transformed_edge_count{0};
    // Ratio of transformed nodes to original nodes (transformed_node_count / original_node_count).
    double node_expansion_ratio{0.0};
    // Path preservation success rate (if applicable, e.g., for ExpansionTransform).
    double path_preservation_rate{0.0};
    // Number of nodes not covered by paths (for PathTraversalTransform).
    std::size_t uncovered_node_count{0};
};

class VGTransform {
public:
    virtual ~VGTransform() = default;

    // Applies the transformation from an ImportedGraph to an AlnGraph.
    // graph_k: The k-mer size used to construct the graph (e.g., in a DBG).
    // pore_k: The k-mer size of the pore model used for squigglization.
    virtual AlnGraph apply(const piru::io::ImportedGraph& imported_graph,
                           std::size_t graph_k,
                           std::size_t pore_k) = 0;

    // Returns statistics about the transformation.
    virtual TransformStats getStats() const { return {}; } // Default empty stats
};

} // namespace piru::index
