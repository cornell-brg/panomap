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
};

// Anchor candidate produced by clustering/chaining.
struct SeedAnchor {
    index::SeedHit target;          // node_id + offset in graph
    std::size_t read_pos{0};        // position in read
    double score{0.0};              // backend-specific score
};

struct ClusterSummary {
    double score{0.0};
    std::vector<SeedAnchor> anchors;  // selected anchors (may be single or multiple)
};

struct SeedClustererConfig {
    std::string backend{"fse"};  // "fse", "probe", "chaining" (future)
    std::size_t max_hash_frequency{0};  // from SeedStore metadata
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
