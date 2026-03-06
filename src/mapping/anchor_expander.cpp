// SPDX-License-Identifier: MIT

#include "mapping/anchor_expander.hpp"

namespace piru::mapping {

// ============================================================================
// PathWalkExpander Implementation
// ============================================================================

PathWalkExpander::PathWalkExpander(const std::vector<std::vector<index::LinearCoordinate>>& coords,
                                   const std::vector<std::size_t>& path_lengths)
    : coords_(coords), path_lengths_(path_lengths) {}

std::vector<PathAnchor> PathWalkExpander::expand(const std::vector<NodeAnchor>& hits) const {
    std::vector<PathAnchor> anchors;

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
            std::int64_t ref_coord = coord.ref_coord + static_cast<std::int64_t>(hit.target.offset);
            std::int64_t anchor_end = ref_coord + static_cast<std::int64_t>(hit.span);

            // Skip anchors that extend past path boundary (cross-node seeds at path end)
            if (coord.path_id < path_lengths_.size()) {
                std::int64_t path_len = static_cast<std::int64_t>(path_lengths_[coord.path_id]);
                if (anchor_end > path_len || ref_coord < 0) {
                    continue;
                }
            }

            PathAnchor anchor;
            anchor.query_pos = hit.read_pos;
            anchor.ref_coord = ref_coord;
            anchor.length = hit.span;
            anchor.path_id = coord.path_id;
            anchor.node_id = node_id;
            anchor.node_offset = hit.target.offset;

            anchors.push_back(anchor);
        }
    }

    return anchors;
}

// ============================================================================
// SuperbubbleExpander Implementation
// ============================================================================

SuperbubbleExpander::SuperbubbleExpander(const index::GraphStore* graph_store)
    : graph_store_(graph_store) {}

std::vector<PathAnchor> SuperbubbleExpander::expand(const std::vector<NodeAnchor>& hits) const {
    std::vector<PathAnchor> anchors;

    // 1:1 mapping for superbubble (no multi-path expansion)
    anchors.reserve(hits.size());

    for (const auto& hit : hits) {
        const std::size_t node_id = hit.target.node_id;

        // Get chain_id and linear_position from GraphStore
        auto chain_id = graph_store_->chainId(node_id);
        auto linear_pos = graph_store_->linearPosition(node_id);

        // Skip nodes without linearization (unmapped/uncovered nodes)
        if (!chain_id.has_value() || !linear_pos.has_value()) {
            continue;
        }

        // Trivial 1:1 mapping to anchor
        PathAnchor anchor;
        anchor.query_pos = hit.read_pos;
        anchor.ref_coord = linear_pos.value() + static_cast<std::int64_t>(hit.target.offset);
        anchor.length = hit.span;
        anchor.path_id = static_cast<std::size_t>(chain_id.value());
        anchor.node_id = node_id;
        anchor.node_offset = hit.target.offset;

        anchors.push_back(anchor);
    }

    return anchors;
}

}  // namespace piru::mapping
