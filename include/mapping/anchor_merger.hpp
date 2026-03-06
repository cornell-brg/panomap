// SPDX-License-Identifier: MIT
// Anchor merging for consolidating adjacent/overlapping anchors.

#pragma once

#include <cstddef>
#include <vector>

#include "mapping/chainer.hpp"

namespace piru::mapping {

// Configuration for anchor merging.
struct AnchorMergerConfig {
    // Currently no parameters needed - merge if overlap, don't merge if gap.
};

// Merge overlapping anchors on the same path.
//
// Algorithm:
// 1. Group anchors by path_id
// 2. Within each path, sort by (ref_coord, query_pos)
// 3. Merge consecutive anchors if they overlap (B starts before A ends)
// 4. Don't merge if there's a gap (B starts after A ends)
// 5. Update length to cover merged span
//
// Benefits:
// - Reduces anchor count for DP chaining
// - More accurate coverage length estimates
// - Better handling of dense seed regions
class AnchorMerger {
public:
    // Merge anchors according to configuration.
    // Returns new vector with merged anchors (input order may not be preserved).
    static std::vector<PathAnchor> merge(const std::vector<PathAnchor>& anchors,
                                     const AnchorMergerConfig& config = AnchorMergerConfig{});
};

}  // namespace piru::mapping
