// SPDX-License-Identifier: MIT
// Seed merging component to consolidate adjacent/overlapping seed hits.

#pragma once

#include <cstddef>
#include <vector>

#include "mapping/chainer.hpp"

namespace piru::mapping {

// Configuration for seed merging behavior.
struct SeedMergerConfig {
    // Maximum gap (in samples) to allow merging between adjacent seeds.
    // Seeds are considered mergeable if:
    // - Same node_id
    // - Gap in query_pos <= merge_tolerance
    // - Gap in node_offset <= merge_tolerance
    std::size_t merge_tolerance{0};
};

// Merges adjacent/overlapping seed hits to reduce redundancy.
//
// Algorithm:
// 1. Sort hits by (node_id, node_offset, query_pos)
// 2. Iterate and merge adjacent hits within tolerance
// 3. Update span field to cover merged region
class SeedMerger {
public:
    // Merge adjacent seed hits within the configured tolerance.
    // Returns a new vector with merged hits (original vector unchanged).
    static std::vector<SeedHitRecord> merge(const std::vector<SeedHitRecord>& hits,
                                            const SeedMergerConfig& config);
};

}  // namespace piru::mapping
