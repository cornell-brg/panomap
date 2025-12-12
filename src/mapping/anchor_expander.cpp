// SPDX-License-Identifier: MIT

#include "mapping/anchor_expander.hpp"

namespace piru::mapping {

AnchorExpander::AnchorExpander(
    const std::vector<std::vector<index::LinearCoordinate>>& coords)
    : coords_(coords) {}

std::vector<Anchor> AnchorExpander::expand(const std::vector<SeedHitRecord>& hits) const {
    std::vector<Anchor> anchors;

    // Pre-allocate approximate space (assume average 1-2 anchors per hit)
    anchors.reserve(hits.size() * 2);

    for (const auto& hit : hits) {
        const std::size_t node_id = hit.target.node_id;

        // Bounds check
        if (node_id >= coords_.size()) {
            continue;  // Skip invalid node IDs
        }

        const auto& node_coords = coords_[node_id];

        // Skip seeds with no linearization coordinates
        if (node_coords.empty()) {
            continue;
        }

        // Emit one anchor per coordinate occurrence
        for (const auto& coord : node_coords) {
            Anchor anchor;
            anchor.query_pos = hit.read_pos;
            anchor.ref_coord = coord.ref_coord + static_cast<std::int64_t>(hit.target.offset);
            anchor.length = hit.span;
            anchor.path_id = coord.path_id;
            anchor.node_id = node_id;
            anchor.node_offset = hit.target.offset;

            anchors.push_back(anchor);
        }
    }

    return anchors;
}

}  // namespace piru::mapping
