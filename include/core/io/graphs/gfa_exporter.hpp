#pragma once

#include <string>

#include "core/index/flat_graph.hpp"
#include "core/io/graphs/graph.hpp"

namespace piru {

class GfaExporter {
public:
  // Export ImportedGraph (bidirectional, base sequences)
  static void dumpImportedGraph(const io::ImportedGraph& graph, const std::string& path);

  // Export FlatGraph (base sequences)
  static void dumpFlatGraph(const index::FlatGraph& graph, const std::string& path);
};

}  // namespace piru
