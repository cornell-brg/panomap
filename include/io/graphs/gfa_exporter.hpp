#pragma once

#include <map>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "io/graphs/graph.hpp"
#include "signal/signal_types.hpp"

namespace piru {

enum class AlnGraphDumpMode { Bases, RawSignal, FuzzyQuantized, AlnQuantized };

class GfaExporter {
public:
    // Export ImportedGraph (bidirectional, base sequences)
    static void dumpImportedGraph(const io::ImportedGraph& graph, const std::string& path);

    // Export AlnGraph with selectable content mode
    // The `data` parameter type depends on the mode:
    // - Bases: nullptr
    // - RawSignal: const std::vector<std::vector<float>>*
    // - FuzzyQuantized: const std::vector<signal::FuzzyQuantizedSignal>*
    // - AlnQuantized: const std::vector<signal::AlignmentQuantizedSignal>*
    static void dumpAlnGraph(const index::AlnGraph& graph, const std::string& path,
                             AlnGraphDumpMode mode = AlnGraphDumpMode::Bases,
                             const void* data = nullptr,
                             const std::map<std::string, std::string>& header_tags = {});
};

}  // namespace piru
