// SPDX-License-Identifier: MIT
// Expands seed hits to anchors using linearization coordinates.

#pragma once

#include <vector>

#include "index/linearizer.hpp"
#include "mapping/anchor.hpp"
#include "mapping/seed_clusterer.hpp"

namespace piru::mapping {

// Expands seed hits to anchors using pre-computed linearization coordinates.
//
// For each seed hit:
// - Lookup linearization coordinates for the node
// - Emit one anchor per coordinate occurrence (handles multiple paths/cycles)
// - Add node_offset to ref_coord to get precise position
//
// Seeds with no linearization coordinates (empty vector) are skipped.
class AnchorExpander {
public:
    // Construct expander with linearization coordinates from a Linearizer.
    // coords[node_id] = vector of LinearCoordinate for that node
    explicit AnchorExpander(const std::vector<std::vector<index::LinearCoordinate>>& coords);

    // Expand seed hits to anchors.
    // Returns vector of anchors (may be larger than input if nodes appear on multiple paths).
    std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const;

private:
    const std::vector<std::vector<index::LinearCoordinate>>& coords_;
};

}  // namespace piru::mapping
