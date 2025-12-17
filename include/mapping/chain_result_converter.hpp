// SPDX-License-Identifier: MIT
// Converts chain results to AlignmentResult for output.

#pragma once

#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "io/results/result.hpp"
#include "mapping/seed_clusterer.hpp"

namespace piru::mapping {

// Configuration for chain-to-result conversion.
struct ChainResultConfig {
    bool primary_only{false};       // Only output primary (best) chain
    std::size_t max_secondary{9};   // Max secondary alignments (default: 9, so up to 10 total)
};

// Converts chain results to AlignmentResult format for PAF/GAF output.
//
// For chain-only mode (no alignment):
// - Query/target coordinates derived from anchor positions
// - MAPQ approximated from chain score
// - No CIGAR (empty mappings or placeholder)
//
// Usage:
//   ChainResultConverter converter(aln_graph, config);
//   auto results = converter.convert(cluster_summary, read_id, read_length);
class ChainResultConverter {
public:
    explicit ChainResultConverter(const index::AlnGraph& graph,
                                  ChainResultConfig config = {});

    // Convert chains to AlignmentResult vector.
    // Returns one result per chain (up to max_secondary + 1).
    std::vector<io::AlignmentResult> convert(
        const ClusterSummary& summary,
        const std::string& read_id,
        std::size_t read_length) const;

private:
    // Convert a single chain (ClusterGroup) to AlignmentResult.
    io::AlignmentResult convertChain(
        const ClusterGroup& chain,
        const std::string& read_id,
        std::size_t read_length,
        bool is_primary) const;

    // Build GAF-style path string from anchors (">node1>node2>node3").
    std::string buildPathString(const std::vector<SeedAnchor>& anchors) const;

    // Get path name from path_id.
    const std::string& getPathName(std::size_t path_id) const;

    // Compute path length (sum of node lengths along path).
    std::size_t computePathLength(std::size_t path_id) const;

    // Estimate MAPQ from chain scores.
    int estimateMapQ(double primary_score, double secondary_score) const;

    const index::AlnGraph& graph_;
    ChainResultConfig config_;

    // Cached path lengths (computed on first access).
    mutable std::vector<std::size_t> path_lengths_;
    mutable bool path_lengths_computed_{false};
};

}  // namespace piru::mapping
