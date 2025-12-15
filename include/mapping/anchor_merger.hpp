// SPDX-License-Identifier: MIT
// Anchor merging for consolidating adjacent/overlapping anchors.

#pragma once

#include <cstddef>
#include <vector>

#include "mapping/anchor.hpp"

namespace piru::mapping {

// Configuration for anchor merging.
struct AnchorMergerConfig {
    // Maximum gap (in samples/bases) to allow merging (default: 0 for exact adjacency).
    std::size_t merge_tolerance{0};
};

// Merge adjacent/overlapping anchors on the same path.
//
// Algorithm:
// 1. Group anchors by path_id
// 2. Within each path, sort by (ref_coord, query_pos)
// 3. Merge consecutive anchors if query gap ≤ tolerance AND ref gap ≤ tolerance
// 4. Update length to cover merged span
//
// Benefits:
// - Reduces anchor count for DP chaining
// - More accurate coverage length estimates
// - Better handling of dense seed regions
class AnchorMerger {
public:
    // Merge anchors according to configuration.
    // Returns new vector with merged anchors (input order may not be preserved).
    static std::vector<Anchor> merge(
        const std::vector<Anchor>& anchors,
        const AnchorMergerConfig& config = AnchorMergerConfig{});
};

}  // namespace piru::mapping
