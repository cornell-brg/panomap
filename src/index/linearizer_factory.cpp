// SPDX-License-Identifier: MIT
// Linearizer factory implementation.

#include "index/linearizer_factory.hpp"

#include <stdexcept>

#include "index/path_walk_linearizer.hpp"
#include "index/superbubble_linearizer.hpp"

namespace piru::index {

std::unique_ptr<Linearizer> make_linearizer(const std::string& backend_name) {
  if (backend_name == "superbubble") {
    return std::make_unique<SuperbubbleLinearizer>();
  } else if (backend_name == "path-walk") {
    return std::make_unique<PathWalkLinearizer>();
  } else {
    throw std::invalid_argument("Unknown linearizer backend: " + backend_name +
                                ". Supported: 'superbubble', 'path-walk'");
  }
}

}  // namespace piru::index
