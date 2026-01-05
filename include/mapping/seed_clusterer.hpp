// SPDX-License-Identifier: MIT
// Interfaces for anchor clustering/chaining backends.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <optional>
#include <vector>

#include "index/seed_store.hpp"

namespace piru::index {
    struct LinearCoordinate;
    class GraphStore;
}

namespace piru::mapping {

struct Anchor;  // Forward declaration (defined in anchor.hpp)

// Minimal hit record used for clustering/chaining.
struct SeedHitRecord {
    index::SeedHit target;      // node_id + offset in graph
    std::size_t read_pos{0};    // seed position in the read
    std::uint64_t hash{0};      // seed hash (for debugging/uniqueness)
    std::size_t span{0};        // coverage length on query (from Seed.length, may be merged)
    std::optional<std::int64_t> chain_id;
    std::optional<std::int64_t> linear_pos;
    std::size_t frequency{0};   // occurrences of this hash in the index
    mutable double score{0.0};  // Computed during clustering (mutable for legacy compatibility)
};

// Anchor candidate produced by clustering/chaining.
struct SeedAnchor {
    index::SeedHit target;          // node_id + offset in graph
    std::size_t read_pos{0};        // position in read
    double score{0.0};              // backend-specific score
    std::size_t cluster_id{0};      // which cluster this anchor belongs to

    // Optional: linear coordinates (for path-walk pipeline debugging)
    std::size_t path_id{0};         // reference path ID
    std::int64_t ref_coord{0};      // linear position on reference path
};

// Group of anchors from a single cluster (for probe-based alignment)
struct ClusterGroup {
    double cluster_score{0.0};          // score of this cluster
    std::vector<SeedAnchor> anchors;    // anchors/probes from this cluster
};

struct ClusterSummary {
    double score{0.0};                      // overall best score
    std::vector<SeedAnchor> anchors;        // flat list (for FSE) - selected/chained anchors
    std::vector<ClusterGroup> clusters;     // grouped by cluster (for Probe)
    std::size_t expanded_anchor_count{0};   // total anchors after expansion (before clustering)
};

struct SeedClustererConfig {
    std::string backend{"fse"};  // "fse", "probe", "dp-chain"
    std::size_t max_hash_frequency{0};  // from SeedStore metadata

    // Clustering parameters (FSE/legacy)
    int diagonal_cutoff{50};           // Max diagonal gap to stay in same cluster
    std::size_t min_cluster_size{2};   // Minimum seeds per cluster
    int max_clusters{-1};              // Max clusters to return (-1 = unlimited)

    // Probe parameters (legacy)
    std::size_t max_probes_per_cluster{10};  // Max probe seeds per cluster
    std::size_t probe_stride{0};             // Override stride (0 = auto-compute from query length)

    // DP chain parameters (used when backend="dp-chain")
    // Tuned for noisy signals (DEV027), original defaults in comments
    std::size_t dp_max_dist{500};               // Max query/ref distance for chaining (banding)
    std::size_t dp_max_diag_dev{500};           // Max diagonal deviation |Δr - Δq|
    double dp_gap_penalty{0.02};                // Penalty per unit gap distance (was: 0.1)
    double dp_diag_penalty{0.05};               // Penalty per unit diagonal deviation (was: 0.5)
    double dp_overlap_penalty{0.90};            // Penalty per unit overlap (was: 2.0)
    double dp_anchor_weight{1.0};               // Weight per anchor length
    std::size_t dp_min_chain_score{0};          // Min score to report a chain (0 = accept any)
    std::size_t dp_max_chains{10};              // Max number of chains to extract
    std::size_t dp_max_skip{25};                // Stop after N consecutive failed chain attempts
    bool dp_merge_chains{true};                 // Merge overlapping chains on same path
};

// Abstract interface for anchor clustering/chaining.
//
// Clusterers operate on anchors (linear reference space) and select
// an optimal subset for alignment extension. All clusterers take the same input
// (std::vector<Anchor>) enabling uniform pipeline architecture.
//
// Different implementations:
// - FSE/Probe: Diagonal clustering (group by path_id, cluster by ref_coord - query_pos)
// - DPChain: Colinear chaining via dynamic programming
// - Noop: Pass-through for debugging
//
// When to use:
// - FSE/Probe: Superbubble pipeline, fast O(n log n) clustering, simple variation graphs
// - DPChain: Path-walk pipeline, complex graphs with cycles, haplotype-aware chaining
//
// Recommended pairings (validated by BatchMapper):
// - Superbubble linearization + FSE/Probe clustering
// - Path-walk linearization + DPChain clustering
//
// Mixing pipelines is allowed but issues a warning (e.g., path-walk + FSE).
class AnchorClusterer {
public:
    virtual ~AnchorClusterer() = default;

    // Cluster/chain anchors and return selected subset.
    // Input: Anchors in linear space (from AnchorExpander)
    // Output: Selected anchors with scores and cluster IDs
    virtual ClusterSummary cluster(const std::vector<Anchor>& anchors) const = 0;

    virtual std::string name() const = 0;
};

using AnchorClustererPtr = std::unique_ptr<AnchorClusterer>;

// Factory function for creating anchor clusterers.
// Note: Expansion is now separate, so linearization_coords is NOT needed here.
// - DiagonalClusterer (FSE/Probe): Operates on anchors directly
// - DPChainClusterer: Operates on anchors directly
// - graph_store parameter kept for potential future use
AnchorClustererPtr make_anchor_clusterer(
    const SeedClustererConfig& config,
    const index::GraphStore* graph_store = nullptr);

// Legacy alias for backward compatibility (deprecated, will be removed in future)
using SeedClusterer = AnchorClusterer;
using SeedClustererPtr = AnchorClustererPtr;
inline AnchorClustererPtr make_seed_clusterer(
    const SeedClustererConfig& config,
    const std::vector<std::vector<index::LinearCoordinate>>* linearization_coords = nullptr,
    const index::GraphStore* graph_store = nullptr) {
    // Ignore linearization_coords - expansion is now separate
    (void)linearization_coords;
    return make_anchor_clusterer(config, graph_store);
}

}  // namespace piru::mapping
