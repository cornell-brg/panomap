#include "io/graphs/gfa_exporter.hpp"

#include <fstream>
#include <iostream>

#include "io/graphs/graph.hpp"

namespace piru {

void GfaExporter::dumpImportedGraph(const io::ImportedGraph& graph, const std::string& path) {
  std::ofstream out(path);
  out << "H\tVN:Z:1.0" << std::endl;

  for (const auto& node : graph.nodes) {
    out << "S\t" << node.id << "\t" << node.sequence << std::endl;
  }

  for (const auto& edge : graph.edges) {
    out << "L\t" << edge.from << "\t" << (edge.from_reverse ? "-" : "+") << "\t" << edge.to << "\t"
        << (edge.to_reverse ? "-" : "+") << "\t";
    if (edge.overlap_bases.has_value()) {
      out << edge.overlap_bases.value() << "M";
    } else {
      out << "*";
    }
    out << std::endl;
  }

  for (const auto& path_entry : graph.paths) {
    out << "P\t" << path_entry.name << "\t";
    for (size_t i = 0; i < path_entry.steps.size(); ++i) {
      if (i > 0) out << ",";
      out << path_entry.steps[i].segment_id << (path_entry.steps[i].is_reverse ? "-" : "+");
    }
    out << "\t*" << std::endl;
  }
}

void GfaExporter::dumpFlatGraph(const index::FlatGraph& graph, const std::string& path) {
  std::ofstream out(path);
  out << "H\tVN:Z:1.0\n";

  for (std::uint32_t i = 0; i < graph.nodeCount(); ++i) {
    out << "S\t" << graph.name(i) << (graph.isReverse(i) ? "-" : "+")
        << "\t" << graph.seq(i) << "\n";
  }

  for (std::uint32_t i = 0; i < graph.nodeCount(); ++i) {
    for (auto it = graph.outBegin(i); it != graph.outEnd(i); ++it) {
      out << "L\t" << i << "\t+\t" << *it << "\t+\t0M\n";
    }
  }

  for (std::uint32_t p = 0; p < graph.pathCount(); ++p) {
    out << "P\t" << graph.pathName(p) << "\t";
    const auto* steps = graph.pathStepsBegin(p);
    std::size_t n = graph.pathStepCount(p);
    for (std::size_t j = 0; j < n; ++j) {
      if (j > 0) out << ",";
      out << steps[j] << "+";
    }
    out << "\t*\n";
  }
}

}  // namespace piru
