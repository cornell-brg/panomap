// SPDX-License-Identifier: MIT
// Factory for creating linearizer backends.

#pragma once

#include <memory>
#include <string>

#include "index/linearizer.hpp"

namespace piru::index {

// Create a linearizer backend by name.
// Supported backends:
// - "superbubble": SuperbubbleLinearizer (uses pseudo-linearization)
// - "path-walk": PathWalkLinearizer (walks reference paths)
//
// Throws std::invalid_argument if backend name is unknown.
std::unique_ptr<Linearizer> make_linearizer(const std::string& backend_name);

}  // namespace piru::index
