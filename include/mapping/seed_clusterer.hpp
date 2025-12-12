// SPDX-License-Identifier: MIT
// Interfaces for seed clustering/chaining backends.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <optional>
#include <vector>

#include "index/seed_store.hpp"

namespace piru::mapping {

// Minimal hit record used for clustering/chaining.
struct SeedHitRecord {
    index::SeedHit target;      // node_id + offset in graph
    std::size_t read_pos{0};    // seed position in the read
    std::uint64_t hash{0};      // seed hash (for debugging/uniqueness)
    std::size_t span{0};        // seed span in read bases (k)
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
};

// Group of anchors from a single cluster (for probe-based alignment)
struct ClusterGroup {
    double cluster_score{0.0};          // score of this cluster
    std::vector<SeedAnchor> anchors;    // anchors/probes from this cluster
};

struct ClusterSummary {
    double score{0.0};                      // overall best score
    std::vector<SeedAnchor> anchors;        // flat list (for FSE)
    std::vector<ClusterGroup> clusters;     // grouped by cluster (for Probe)
};

struct SeedClustererConfig {
    std::string backend{"fse"};  // "fse", "probe", "chaining" (future)
    std::size_t max_hash_frequency{0};  // from SeedStore metadata

    // Clustering parameters (FSE/legacy)
    int diagonal_cutoff{50};           // Max diagonal gap to stay in same cluster
    std::size_t min_cluster_size{2};   // Minimum seeds per cluster
    int max_clusters{-1};              // Max clusters to return (-1 = unlimited)

    // Probe parameters (legacy)
    std::size_t max_probes_per_cluster{10};  // Max probe seeds per cluster
    std::size_t probe_stride{0};             // Override stride (0 = auto-compute from query length)
};

class SeedClusterer {
public:
    virtual ~SeedClusterer() = default;

    // Process seed hits for one read and return selected anchors/summary.
    virtual ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const = 0;
    virtual std::string name() const = 0;
};

using SeedClustererPtr = std::unique_ptr<SeedClusterer>;

SeedClustererPtr make_seed_clusterer(const SeedClustererConfig& config);

}  // namespace piru::mapping
