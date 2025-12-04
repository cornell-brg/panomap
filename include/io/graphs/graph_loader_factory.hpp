// SPDX-License-Identifier: MIT
// Factory helpers to select a graph loader based on format/path.

#pragma once

#include <memory>
#include <string>

#include "io/graphs/graph_loader.hpp"

namespace piru::io {

// Attempt to choose a loader based on path/extension.
// Returns nullptr if no suitable loader is available.
GraphLoaderPtr make_graph_loader(const std::string& path);

}  // namespace piru::io
