// SPDX-License-Identifier: MIT
// Superbubble-based linearization using pseudo-linearization algorithm.

#pragma once

#include "index/linearizer.hpp"

namespace piru::index {

// Linearizer backend using superbubble-based pseudo-linearization.
// Assigns chain IDs and local linear coordinates to nodes within superbubble chains.
// Nodes not assigned to any chain (complex cyclic regions, tips) get no coordinates.
class SuperbubbleLinearizer : public Linearizer {
public:
  SuperbubbleLinearizer() = default;

  std::vector<std::vector<LinearCoordinate>> linearize(const AlnGraph& graph) const override;

  std::string name() const override { return "superbubble"; }
};

}  // namespace piru::index
