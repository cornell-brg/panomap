// SPDX-License-Identifier: MIT
// Interface for loading graph representations from files.

#pragma once

#include <memory>
#include <string>

#include "io/graphs/graph.hpp"

namespace piru::io {

class GraphLoader {
public:
  virtual ~GraphLoader() = default;

  // Populate `graph` from the source file. Returns false on parsing/open errors.
  virtual bool load(ImportedGraph& graph) = 0;

  // Human-readable format name (e.g., gfa, vg).
  virtual std::string get_format_name() const = 0;
};

using GraphLoaderPtr = std::unique_ptr<GraphLoader>;

}  // namespace piru::io
