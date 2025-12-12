// SPDX-License-Identifier: MIT
// Path-guided linearization by walking reference paths.

#pragma once

#include "index/linearizer.hpp"

namespace piru::index {

// Linearizer backend that walks reference paths to assign global linear coordinates.
// Each node gets multiple coordinates (one per path occurrence).
// Handles cycles: if a node appears multiple times on same path, each occurrence
// gets a separate coordinate.
// Throws std::runtime_error if graph has no paths.
class PathWalkLinearizer : public Linearizer {
public:
  PathWalkLinearizer() = default;

  std::vector<std::vector<LinearCoordinate>> linearize(const AlnGraph& graph) const override;

  std::string name() const override { return "path-walk"; }
};

}  // namespace piru::index
