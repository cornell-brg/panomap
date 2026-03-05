// SPDX-License-Identifier: MIT
// Linear coordinate type for graph node linearization.

#pragma once

#include <cstddef>
#include <cstdint>

namespace piru::index {

// Linear coordinate assigned to a node.
// A node may have multiple coordinates if it appears on multiple paths
// or appears multiple times on the same path (cycles).
struct LinearCoordinate {
    std::size_t path_id;     // Which reference path
    std::int64_t ref_coord;  // Linear position along that path

    LinearCoordinate(std::size_t pid, std::int64_t coord) : path_id(pid), ref_coord(coord) {}
};

}  // namespace piru::index
