// SPDX-License-Identifier: MIT
#include "index/path_guided_transform.hpp"
#include <algorithm>
#include <stdexcept>

namespace piru::index {

//------------------------------------------------------------------------------
// Task 1: Transform ImportedGraph to AlnGraph (node splitting)
//------------------------------------------------------------------------------
std::pair<AlnGraph, NodeMapping>
PathGuidedTransform::importedGraphToAlnGraph(const piru::io::ImportedGraph& imported) {
  AlnGraph graph;
  NodeMapping node_mapping;

  // Step 1: Create forward and reverse variants for each node
  for (const auto& imported_node : imported.nodes) {
    // Forward variant
    AlnNode fwd_node;
    fwd_node.label = imported_node.id + "_F";
    fwd_node.original_id = imported_node.id;
    fwd_node.is_reverse = false;
    fwd_node.sequence = imported_node.sequence;  // No context yet
    std::size_t fwd_id = graph.addNode(fwd_node);

    // Reverse variant
    AlnNode rev_node;
    rev_node.label = imported_node.id + "_R";
    rev_node.original_id = imported_node.id;
    rev_node.is_reverse = true;
    rev_node.sequence = revcomp(imported_node.sequence);  // No context yet
    std::size_t rev_id = graph.addNode(rev_node);

    // Store mapping
    node_mapping[imported_node.id] = {fwd_id, rev_id};
  }

  // Step 2: Transform edges to connect directional nodes
  for (const auto& imported_edge : imported.edges) {
    const auto from_it = node_mapping.find(imported_edge.from);
    const auto to_it = node_mapping.find(imported_edge.to);

    if (from_it == node_mapping.end() || to_it == node_mapping.end()) {
      continue;  // Skip edges referencing missing nodes
    }

    // Determine which directional nodes to connect based on reverse flags
    // ImportedGraphEdge has: from_reverse and to_reverse
    // If from_reverse = false → use forward variant of 'from'
    // If from_reverse = true → use reverse variant of 'from'
    // Same logic for 'to'

    std::size_t from_aln_id = imported_edge.from_reverse
                                  ? from_it->second.second  // reverse variant
                                  : from_it->second.first;  // forward variant

    std::size_t to_aln_id = imported_edge.to_reverse
                                ? to_it->second.second  // reverse variant
                                : to_it->second.first;  // forward variant

    AlnEdge edge;
    edge.from = from_aln_id;
    edge.to = to_aln_id;
    edge.overlap_bases = 0;  // No overlap yet (context will be added in Task 2)

    graph.addEdge(edge);
  }

  return {graph, node_mapping};
}

//------------------------------------------------------------------------------
// Task 2 & 3: Walk paths and add k-1 context (includes coverage tracking)
//------------------------------------------------------------------------------
void PathGuidedTransform::walkPathsAndAddContext(
    AlnGraph& graph,
    const piru::io::ImportedGraph& imported,
    const NodeMapping& node_mapping,
    std::size_t pore_k,
    std::unordered_set<std::string>& covered_nodes) {

  const std::size_t k_minus_1 = pore_k - 1;

  // Track node variants: (original_id, is_reverse, context) → aln_node_id
  std::unordered_map<std::string, std::size_t> variant_map;

  for (const auto& path : imported.paths) {
    AlnPath aln_path;
    aln_path.name = path.name;

    for (std::size_t i = 0; i < path.steps.size(); ++i) {
      const auto& step = path.steps[i];
      const std::string& node_id = step.segment_id;
      const bool is_reverse = step.is_reverse;

      // Mark this node as covered
      covered_nodes.insert(node_id);

      // Get the base AlnGraph node (before adding context)
      const auto mapping_it = node_mapping.find(node_id);
      if (mapping_it == node_mapping.end()) {
        continue;  // Skip if node not found
      }

      std::size_t base_aln_id = is_reverse ? mapping_it->second.second
                                           : mapping_it->second.first;

      // Determine successor context
      std::string context;
      std::size_t successor_aln_id = 0;

      if (i + 1 < path.steps.size()) {
        // Get successor from next step
        const auto& succ_step = path.steps[i + 1];
        const auto succ_mapping_it = node_mapping.find(succ_step.segment_id);

        if (succ_mapping_it != node_mapping.end()) {
          successor_aln_id = succ_step.is_reverse ? succ_mapping_it->second.second
                                                   : succ_mapping_it->second.first;

          const AlnNode& succ_node = graph.node(successor_aln_id);
          context = getKMinus1Context(succ_node.sequence, k_minus_1);
        }
      }
      // If no successor (path end), context remains empty

      // Create variant key: (original_id, is_reverse, context)
      std::string variant_key = node_id + "|" +
                                (is_reverse ? "R" : "F") + "|" +
                                context;

      std::size_t node_with_context_id;

      auto variant_it = variant_map.find(variant_key);
      if (variant_it != variant_map.end()) {
        // Variant already exists, reuse it
        node_with_context_id = variant_it->second;
      } else {
        // Create new variant with context
        const AlnNode& base_node = graph.node(base_aln_id);

        AlnNode node_with_context;
        node_with_context.label = base_node.label + "_ctx";
        node_with_context.original_id = base_node.original_id;
        node_with_context.is_reverse = base_node.is_reverse;
        node_with_context.sequence = base_node.sequence + context;

        node_with_context_id = graph.addNode(node_with_context);
        variant_map[variant_key] = node_with_context_id;

        // Wire edges: this variant should point to its successor
        if (successor_aln_id != 0) {
          AlnEdge edge;
          edge.from = node_with_context_id;
          edge.to = successor_aln_id;
          edge.overlap_bases = k_minus_1;
          graph.addEdge(edge);
        }
      }

      // Add step to path
      AlnPathStep aln_step;
      aln_step.node_id = std::to_string(node_with_context_id);
      aln_step.is_reverse = is_reverse;
      aln_path.steps.push_back(aln_step);
    }

    // Add reconstructed path to graph
    graph.addPath(aln_path);
  }
}

//------------------------------------------------------------------------------
// Task 4: Orchestrate in apply()
//------------------------------------------------------------------------------
AlnGraph PathGuidedTransform::apply(const piru::io::ImportedGraph& imported,
                                     std::size_t graph_k,
                                     std::size_t pore_k) {
  // Task 1: Transform to directional AlnGraph
  auto [graph, node_mapping] = importedGraphToAlnGraph(imported);

  // Task 2 & 3: Walk paths, add context, track coverage
  std::unordered_set<std::string> covered_nodes;
  walkPathsAndAddContext(graph, imported, node_mapping, pore_k, covered_nodes);

  // Compute statistics
  stats_.original_node_count = imported.nodes.size();
  stats_.transformed_node_count = graph.nodeCount();
  stats_.original_edge_count = imported.edges.size();
  stats_.transformed_edge_count = graph.edgeCount();
  stats_.node_expansion_ratio = stats_.original_node_count > 0
                                    ? static_cast<double>(stats_.transformed_node_count) /
                                          static_cast<double>(stats_.original_node_count)
                                    : 0.0;
  stats_.uncovered_node_count = imported.nodes.size() - covered_nodes.size();

  // TODO Stage 3: Handle uncovered nodes with expansion
  // For now, uncovered nodes are simply not included in the output

  return graph;
}

//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
std::string PathGuidedTransform::revcomp(const std::string& seq) const {
  static const std::unordered_map<char, char> comp = {
      {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
      {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
      {'N', 'N'}, {'n', 'n'}};

  std::string rc;
  rc.reserve(seq.size());
  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    auto comp_it = comp.find(*it);
    rc += (comp_it != comp.end()) ? comp_it->second : 'N';
  }
  return rc;
}

std::string PathGuidedTransform::getKMinus1Context(const std::string& seq,
                                                    std::size_t k_minus_1) const {
  if (seq.size() < k_minus_1) {
    return seq;  // Return full sequence if shorter than k-1
  }
  return seq.substr(0, k_minus_1);
}

} // namespace piru::index
