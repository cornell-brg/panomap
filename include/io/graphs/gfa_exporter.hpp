#pragma once

#include <string>

#include "index/aln_graph.hpp"
#include "io/graphs/graph.hpp"

namespace piru {

enum class AlnGraphDumpMode {
  Bases,
  RawSignal,
  FuzzyQuantized,
  AlnQuantized
};

class GfaExporter {
 public:
  // Export ImportedGraph (bidirectional, base sequences)
  static void dumpImportedGraph(const io::ImportedGraph& graph,
                                const std::string& path);

  // Export AlnGraph with selectable content mode
  static void dumpAlnGraph(const index::AlnGraph& graph,
                           const std::string& path,
                           AlnGraphDumpMode mode = AlnGraphDumpMode::Bases,
                           const void* data = nullptr);
};

}  // namespace piru
