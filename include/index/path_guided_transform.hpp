// SPDX-License-Identifier: MIT
// Path-guided VG transformation: uses embedded haplotype paths to add k-1 context.

#pragma once

#include "io/graphs/graph.hpp"
#include "index/aln_graph.hpp"
#include "index/vg_transform.hpp"
#include <unordered_map>
#include <unordered_set>
#include <string>

namespace piru::index {

// Mapping from original node ID to (forward_aln_id, reverse_aln_id)
using NodeMapping = std::unordered_map<std::string, std::pair<std::size_t, std::size_t>>;

// Structure to track context information for a node
struct ContextInfo {
  std::string context;           // k-1 bases from successor
  std::size_t successor_aln_id;  // Which AlnGraph node provided this context
};

class PathGuidedTransform : public VGTransform {
public:
  PathGuidedTransform() = default;
  ~PathGuidedTransform() override = default;

  AlnGraph apply(const piru::io::ImportedGraph& imported,
                 std::size_t graph_k,
                 std::size_t pore_k) override;

  TransformStats getStats() const override { return stats_; }

private:
  // Task 1: Transform ImportedGraph to AlnGraph (node splitting, no context yet)
  std::pair<AlnGraph, NodeMapping>
  importedGraphToAlnGraph(const piru::io::ImportedGraph& imported);

  // Task 2: Walk paths and add k-1 context from successors
  void walkPathsAndAddContext(AlnGraph& graph,
                               const piru::io::ImportedGraph& imported,
                               const NodeMapping& node_mapping,
                               std::size_t pore_k,
                               std::unordered_set<std::string>& covered_nodes);

  // Stage 3: Expand uncovered nodes by traversing successors
  void expandUncoveredNodes(AlnGraph& graph,
                            const std::unordered_set<std::string>& uncovered_node_ids,
                            const NodeMapping& node_mapping,
                            std::size_t k_minus_1);

  // Helper: Collect all k-1 contexts via depth-limited traversal
  std::vector<ContextInfo> collectKMinus1Contexts(const AlnGraph& graph,
                                                   std::size_t start_node_id,
                                                   std::size_t depth,
                                                   std::unordered_set<std::size_t>& visited) const;

  // Helper: Get reverse complement of a sequence
  std::string revcomp(const std::string& seq) const;

  // Helper: Get k-1 bases from a node's sequence
  std::string getKMinus1Context(const std::string& seq, std::size_t k_minus_1) const;

  TransformStats stats_;
};

} // namespace piru::index
