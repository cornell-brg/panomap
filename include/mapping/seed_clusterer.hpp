// SPDX-License-Identifier: MIT
// Interfaces for anchor clustering/chaining backends.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "index/seed_store.hpp"

namespace piru::mapping {

struct Anchor;  // Forward declaration (defined in anchor.hpp)

// Minimal hit record used for clustering/chaining.
struct SeedHitRecord {
    index::SeedHit target;    // node_id + offset in graph
    std::size_t read_pos{0};  // seed position in the read
    std::uint64_t hash{0};    // seed hash (for debugging/uniqueness)
    std::size_t span{0};      // coverage length on query (from Seed.length, may be merged)
    std::optional<std::int64_t> chain_id;
    std::optional<std::int64_t> linear_pos;
    std::size_t frequency{0};   // occurrences of this hash in the index
    mutable double score{0.0};  // Computed during clustering (mutable for legacy compatibility)
};

// Anchor candidate produced by clustering/chaining.
struct SeedAnchor {
    index::SeedHit target;      // node_id + offset in graph
    std::size_t read_pos{0};    // position in read
    double score{0.0};          // backend-specific score
    std::size_t cluster_id{0};  // which cluster this anchor belongs to

    // Optional: linear coordinates (for path-walk pipeline debugging)
    std::size_t path_id{0};     // reference path ID
    std::int64_t ref_coord{0};  // linear position on reference path
};

// Group of anchors from a single cluster.
struct ClusterGroup {
    double cluster_score{0.0};
    std::vector<SeedAnchor> anchors;
};

struct ClusterSummary {
    double score{0.0};
    std::vector<SeedAnchor> anchors;       // flat list of selected/chained anchors
    std::vector<ClusterGroup> clusters;    // grouped by cluster
    std::size_t expanded_anchor_count{0};  // total anchors before clustering
};

struct SeedClustererConfig {
    std::string backend{"dp-chain"};
    std::size_t max_hash_frequency{0};  // from SeedStore metadata

    // DP chain parameters (tuned for noisy signals, DEV027)
    std::size_t dp_max_dist{500};        // Max query/ref distance for chaining (banding)
    std::size_t dp_max_diag_dev{500};    // Max diagonal deviation |dr - dq|
    double dp_gap_penalty{0.02};         // Penalty per unit gap distance
    double dp_diag_penalty{0.05};        // Penalty per unit diagonal deviation
    double dp_overlap_penalty{0.90};     // Penalty per unit overlap
    double dp_anchor_weight{1.0};        // Weight per anchor length
    std::size_t dp_min_chain_score{12};  // Min score to report a chain
    std::size_t dp_max_chains{10};       // Max number of chains to extract
    std::size_t dp_max_skip{25};         // Stop after N consecutive failed chain attempts
    bool dp_merge_chains{false};         // Merge overlapping chains on same path
};

// Abstract interface for anchor clustering/chaining.
// Clusterers operate on anchors (linear reference space) and select
// an optimal subset for alignment extension.
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
AnchorClustererPtr make_anchor_clusterer(const SeedClustererConfig& config);

}  // namespace piru::mapping
