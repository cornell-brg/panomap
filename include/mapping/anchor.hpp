// SPDX-License-Identifier: MIT
// Anchor data structure for DP-based seed chaining.

#pragma once

#include <cstddef>
#include <cstdint>

namespace piru::mapping {

// Anchor represents a seed hit mapped to a linear reference coordinate.
//
// Anchors are produced by expanding seed hits using linearization coordinates.
// One seed hit may expand to multiple anchors (one per path occurrence of the node).
struct PathAnchor {
    std::size_t query_pos{0};   // Position in query/read
    std::int64_t ref_coord{0};  // Linear position along reference path
    std::size_t length{0};      // Coverage length (from seed span)
    std::size_t path_id{0};     // Which reference path this anchor belongs to

    // Back-references for debugging and alignment construction
    std::size_t node_id{0};      // Graph node this anchor came from
    std::size_t node_offset{0};  // Offset within the node
};

}  // namespace piru::mapping
