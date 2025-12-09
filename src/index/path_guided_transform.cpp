// SPDX-License-Identifier: MIT
#include "index/path_guided_transform.hpp"
#include <algorithm>
#include <stdexcept>

namespace piru::index {

//------------------------------------------------------------------------------
// Task 1: Setup node mapping (nodes created later with context in Task 2)
//------------------------------------------------------------------------------
std::pair<AlnGraph, NodeMapping>
PathGuidedTransform::importedGraphToAlnGraph(const piru::io::ImportedGraph& imported) {
  AlnGraph graph;
  NodeMapping node_mapping;

  // Don't create nodes yet - they will be created with context in Task 2
  // Just initialize node_mapping with placeholders (will be populated in Task 2)
  for (const auto& imported_node : imported.nodes) {
    node_mapping[imported_node.id] = {0, 0};  // Placeholder: will be set in Task 2
  }

  return {graph, node_mapping};
}

//------------------------------------------------------------------------------
// Task 2 & 3: Walk paths and add k-1 context (includes coverage tracking)
//------------------------------------------------------------------------------
std::unordered_map<std::string, std::size_t> PathGuidedTransform::walkPathsAndAddContext(
    AlnGraph& graph,
    const piru::io::ImportedGraph& imported,
    NodeMapping& node_mapping,  // Non-const: we update it as we create nodes
    std::size_t pore_k,
    std::unordered_set<std::string>& covered_nodes) {

  const std::size_t k_minus_1 = pore_k - 1;

  // Track node variants: (original_id, is_reverse, context) → aln_node_id
  std::unordered_map<std::string, std::size_t> variant_map;

  // Get imported node by id
  std::unordered_map<std::string, const io::ImportedGraphNode*> imported_node_map;
  for (const auto& node : imported.nodes) {
    imported_node_map[node.id] = &node;
  }

  for (const auto& path : imported.paths) {
    AlnPath aln_path;
    aln_path.name = path.name;

    for (std::size_t i = 0; i < path.steps.size(); ++i) {
      const auto& step = path.steps[i];
      const std::string& node_id = step.segment_id;
      const bool is_reverse = step.is_reverse;

      // Mark this node as covered
      covered_nodes.insert(node_id);

      // Get the imported node
      const auto imported_it = imported_node_map.find(node_id);
      if (imported_it == imported_node_map.end()) {
        continue;
      }
      const io::ImportedGraphNode& imported_node = *imported_it->second;

      // Determine successor context
      std::string context;

      if (i + 1 < path.steps.size()) {
        // Get successor from next step
        const auto& succ_step = path.steps[i + 1];
        const auto succ_imported_it = imported_node_map.find(succ_step.segment_id);

        if (succ_imported_it != imported_node_map.end()) {
          const io::ImportedGraphNode& succ_imported = *succ_imported_it->second;
          // Context comes from the beginning of successor sequence
          // Handle orientation: if successor is reverse, take context from its revcomp
          std::string succ_seq = succ_step.is_reverse ? revcomp(succ_imported.sequence)
                                                       : succ_imported.sequence;
          context = getKMinus1Context(succ_seq, k_minus_1);
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
        // Create new node with context
        AlnNode node;
        node.label = node_id + (is_reverse ? "_R" : "_F");
        node.original_id = node_id;
        node.is_reverse = is_reverse;
        // Apply orientation to sequence and append context
        std::string base_seq = is_reverse ? revcomp(imported_node.sequence)
                                           : imported_node.sequence;
        node.sequence = base_seq + context;

        node_with_context_id = graph.addNode(node);
        variant_map[variant_key] = node_with_context_id;

        // Update node_mapping for first occurrence of F/R variant
        auto& mapping = node_mapping[node_id];
        if (is_reverse) {
          if (mapping.second == 0) mapping.second = node_with_context_id;
        } else {
          if (mapping.first == 0) mapping.first = node_with_context_id;
        }
      }

      // Add step to path
      AlnPathStep aln_step;
      aln_step.node_id = graph.node(node_with_context_id).label;
      aln_step.is_reverse = is_reverse;
      aln_path.steps.push_back(aln_step);
    }

    // Add reconstructed path to graph
    graph.addPath(aln_path);
  }

  // Walk paths in reverse direction to populate reverse nodes
  for (const auto& path : imported.paths) {
    AlnPath reverse_aln_path;
    reverse_aln_path.name = path.name + "_reverse";

    // Iterate through path steps in reverse order with flipped orientations
    for (std::size_t i = path.steps.size(); i-- > 0; ) {
      const auto& step = path.steps[i];
      const std::string& node_id = step.segment_id;
      const bool is_reverse = !step.is_reverse;  // Flip orientation

      // Mark this node as covered
      covered_nodes.insert(node_id);

      // Get the imported node
      const auto imported_it = imported_node_map.find(node_id);
      if (imported_it == imported_node_map.end()) {
        continue;
      }
      const io::ImportedGraphNode& imported_node = *imported_it->second;

      // Determine successor context (from next step in reverse walk)
      std::string context;

      if (i > 0) {
        // Get predecessor in original path (successor in reverse path)
        const auto& succ_step = path.steps[i - 1];
        const auto succ_imported_it = imported_node_map.find(succ_step.segment_id);

        if (succ_imported_it != imported_node_map.end()) {
          const io::ImportedGraphNode& succ_imported = *succ_imported_it->second;
          // Flip orientation for successor as well
          bool succ_is_reverse = !succ_step.is_reverse;
          std::string succ_seq = succ_is_reverse ? revcomp(succ_imported.sequence)
                                                   : succ_imported.sequence;
          context = getKMinus1Context(succ_seq, k_minus_1);
        }
      }

      // Create variant key
      std::string variant_key = node_id + "|" +
                                (is_reverse ? "R" : "F") + "|" +
                                context;

      std::size_t node_with_context_id;

      auto variant_it = variant_map.find(variant_key);
      if (variant_it != variant_map.end()) {
        // Variant already exists, reuse it
        node_with_context_id = variant_it->second;
      } else {
        // Create new node with context
        AlnNode node;
        node.label = node_id + (is_reverse ? "_R" : "_F");
        node.original_id = node_id;
        node.is_reverse = is_reverse;
        std::string base_seq = is_reverse ? revcomp(imported_node.sequence)
                                           : imported_node.sequence;
        node.sequence = base_seq + context;

        node_with_context_id = graph.addNode(node);
        variant_map[variant_key] = node_with_context_id;

        // Update node_mapping for first occurrence
        auto& mapping = node_mapping[node_id];
        if (is_reverse) {
          if (mapping.second == 0) mapping.second = node_with_context_id;
        } else {
          if (mapping.first == 0) mapping.first = node_with_context_id;
        }
      }

      // Add step to reverse path
      AlnPathStep aln_step;
      aln_step.node_id = graph.node(node_with_context_id).label;
      aln_step.is_reverse = is_reverse;
      reverse_aln_path.steps.push_back(aln_step);
    }

    // Add reverse path to graph
    graph.addPath(reverse_aln_path);
  }

  // Now create edges between nodes based on path order and context
  for (const auto& path : imported.paths) {
    for (std::size_t i = 0; i + 1 < path.steps.size(); ++i) {
      const auto& step = path.steps[i];
      const auto& succ_step = path.steps[i + 1];

      // Get context from successor
      const auto succ_imported_it = imported_node_map.find(succ_step.segment_id);
      if (succ_imported_it == imported_node_map.end()) continue;

      std::string succ_seq = succ_step.is_reverse
                                 ? revcomp(succ_imported_it->second->sequence)
                                 : succ_imported_it->second->sequence;
      std::string context = getKMinus1Context(succ_seq, k_minus_1);

      // Build variant keys
      std::string from_key = step.segment_id + "|" +
                             (step.is_reverse ? "R" : "F") + "|" + context;
      std::string to_key = succ_step.segment_id + "|" +
                           (succ_step.is_reverse ? "R" : "F") + "|" + "";

      auto from_it = variant_map.find(from_key);
      auto to_it = variant_map.find(to_key);

      if (from_it != variant_map.end() && to_it != variant_map.end()) {
        // Check if edge already exists
        bool edge_exists = false;
        for (const auto& existing_edge : graph.edges()) {
          if (existing_edge.from == from_it->second && existing_edge.to == to_it->second) {
            edge_exists = true;
            break;
          }
        }

        if (!edge_exists) {
          AlnEdge edge;
          edge.from = from_it->second;
          edge.to = to_it->second;
          edge.overlap_bases = k_minus_1;
          graph.addEdge(edge);
        }
      }
    }
  }

  // Create edges for reverse paths
  for (const auto& path : imported.paths) {
    for (std::size_t i = path.steps.size() - 1; i > 0; --i) {
      const auto& step = path.steps[i];
      const auto& succ_step = path.steps[i - 1];

      // Flip orientations for reverse walk
      bool step_is_reverse = !step.is_reverse;
      bool succ_is_reverse = !succ_step.is_reverse;

      // Get context from successor
      const auto succ_imported_it = imported_node_map.find(succ_step.segment_id);
      if (succ_imported_it == imported_node_map.end()) continue;

      std::string succ_seq = succ_is_reverse
                                 ? revcomp(succ_imported_it->second->sequence)
                                 : succ_imported_it->second->sequence;
      std::string context = getKMinus1Context(succ_seq, k_minus_1);

      // Build variant keys
      std::string from_key = step.segment_id + "|" +
                             (step_is_reverse ? "R" : "F") + "|" + context;
      std::string to_key = succ_step.segment_id + "|" +
                           (succ_is_reverse ? "R" : "F") + "|" + "";

      auto from_it = variant_map.find(from_key);
      auto to_it = variant_map.find(to_key);

      if (from_it != variant_map.end() && to_it != variant_map.end()) {
        // Check if edge already exists
        bool edge_exists = false;
        for (const auto& existing_edge : graph.edges()) {
          if (existing_edge.from == from_it->second && existing_edge.to == to_it->second) {
            edge_exists = true;
            break;
          }
        }

        if (!edge_exists) {
          AlnEdge edge;
          edge.from = from_it->second;
          edge.to = to_it->second;
          edge.overlap_bases = k_minus_1;
          graph.addEdge(edge);
        }
      }
    }
  }

  return variant_map;
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
  auto variant_map = walkPathsAndAddContext(graph, imported, node_mapping, pore_k, covered_nodes);

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

  // Stage 3: Handle uncovered nodes with expansion
  std::unordered_set<std::string> uncovered_node_ids;
  for (const auto& node : imported.nodes) {
    if (covered_nodes.find(node.id) == covered_nodes.end()) {
      uncovered_node_ids.insert(node.id);
    }
  }

  if (!uncovered_node_ids.empty()) {
    const std::size_t k_minus_1 = pore_k - 1;
    expandUncoveredNodes(graph, uncovered_node_ids, node_mapping, k_minus_1, imported);
  }

  // Stage 4: Create edges for all imported graph connections
  // Use variant_map to find nodes with proper context, fall back to base nodes
  const std::size_t k_minus_1 = pore_k - 1;

  // Build imported node map for context lookup
  std::unordered_map<std::string, const io::ImportedGraphNode*> imported_node_map;
  for (const auto& node : imported.nodes) {
    imported_node_map[node.id] = &node;
  }

  for (const auto& imported_edge : imported.edges) {
    const auto from_imported_it = imported_node_map.find(imported_edge.from);
    const auto to_imported_it = imported_node_map.find(imported_edge.to);

    if (from_imported_it == imported_node_map.end() || to_imported_it == imported_node_map.end()) {
      continue;
    }

    const io::ImportedGraphNode& from_imported = *from_imported_it->second;
    const io::ImportedGraphNode& to_imported = *to_imported_it->second;

    // Determine orientations
    bool from_is_reverse = imported_edge.from_reverse;
    bool to_is_reverse = imported_edge.to_reverse;

    // Get context from successor
    std::string to_seq = to_is_reverse ? revcomp(to_imported.sequence) : to_imported.sequence;
    std::string context = getKMinus1Context(to_seq, k_minus_1);

    // Build variant key for from node (with context to successor)
    std::string from_key = imported_edge.from + "|" +
                           (from_is_reverse ? "R" : "F") + "|" + context;

    // Build variant key for to node (no context - it's the successor)
    std::string to_key = imported_edge.to + "|" +
                         (to_is_reverse ? "R" : "F") + "|" + "";

    // Try to find variants in variant_map
    auto from_variant_it = variant_map.find(from_key);
    auto to_variant_it = variant_map.find(to_key);

    // Fall back to node_mapping if not in variant_map
    std::size_t from_aln_id = 0;
    std::size_t to_aln_id = 0;

    if (from_variant_it != variant_map.end()) {
      from_aln_id = from_variant_it->second;
    } else {
      // Use base node from node_mapping
      auto from_mapping_it = node_mapping.find(imported_edge.from);
      if (from_mapping_it != node_mapping.end()) {
        from_aln_id = from_is_reverse ? from_mapping_it->second.second
                                       : from_mapping_it->second.first;
      }
    }

    if (to_variant_it != variant_map.end()) {
      to_aln_id = to_variant_it->second;
    } else {
      // Use base node from node_mapping
      auto to_mapping_it = node_mapping.find(imported_edge.to);
      if (to_mapping_it != node_mapping.end()) {
        to_aln_id = to_is_reverse ? to_mapping_it->second.second
                                   : to_mapping_it->second.first;
      }
    }

    if (from_aln_id == 0 || to_aln_id == 0) {
      continue;  // Skip if nodes don't exist
    }

    // Check if edge already exists
    bool edge_exists = false;
    for (const auto& existing_edge : graph.edges()) {
      if (existing_edge.from == from_aln_id && existing_edge.to == to_aln_id) {
        edge_exists = true;
        break;
      }
    }

    if (!edge_exists) {
      AlnEdge edge;
      edge.from = from_aln_id;
      edge.to = to_aln_id;
      edge.overlap_bases = k_minus_1;
      graph.addEdge(edge);
    }

    // Reverse complement edge: to_reverse → from_reverse
    bool rev_from_is_reverse = !to_is_reverse;
    bool rev_to_is_reverse = !from_is_reverse;

    // Context for reverse edge (from what's now the predecessor)
    std::string rev_to_seq = rev_to_is_reverse ? revcomp(from_imported.sequence) : from_imported.sequence;
    std::string rev_context = getKMinus1Context(rev_to_seq, k_minus_1);

    std::string rev_from_key = imported_edge.to + "|" +
                               (rev_from_is_reverse ? "R" : "F") + "|" + rev_context;
    std::string rev_to_key = imported_edge.from + "|" +
                             (rev_to_is_reverse ? "R" : "F") + "|" + "";

    auto rev_from_variant_it = variant_map.find(rev_from_key);
    auto rev_to_variant_it = variant_map.find(rev_to_key);

    std::size_t rev_from_aln_id = 0;
    std::size_t rev_to_aln_id = 0;

    if (rev_from_variant_it != variant_map.end()) {
      rev_from_aln_id = rev_from_variant_it->second;
    } else {
      auto rev_from_mapping_it = node_mapping.find(imported_edge.to);
      if (rev_from_mapping_it != node_mapping.end()) {
        rev_from_aln_id = rev_from_is_reverse ? rev_from_mapping_it->second.second
                                               : rev_from_mapping_it->second.first;
      }
    }

    if (rev_to_variant_it != variant_map.end()) {
      rev_to_aln_id = rev_to_variant_it->second;
    } else {
      auto rev_to_mapping_it = node_mapping.find(imported_edge.from);
      if (rev_to_mapping_it != node_mapping.end()) {
        rev_to_aln_id = rev_to_is_reverse ? rev_to_mapping_it->second.second
                                           : rev_to_mapping_it->second.first;
      }
    }

    if (rev_from_aln_id == 0 || rev_to_aln_id == 0) {
      continue;
    }

    edge_exists = false;
    for (const auto& existing_edge : graph.edges()) {
      if (existing_edge.from == rev_from_aln_id && existing_edge.to == rev_to_aln_id) {
        edge_exists = true;
        break;
      }
    }

    if (!edge_exists) {
      AlnEdge rev_edge;
      rev_edge.from = rev_from_aln_id;
      rev_edge.to = rev_to_aln_id;
      rev_edge.overlap_bases = k_minus_1;
      graph.addEdge(rev_edge);
    }
  }

  return graph;
}

//------------------------------------------------------------------------------
// Stage 3: Expansion algorithm for uncovered nodes
//------------------------------------------------------------------------------

// Recursive depth-limited traversal to collect all k-1 contexts
// Returns contexts by greedily collecting bases from successor paths
std::vector<ContextInfo> PathGuidedTransform::collectKMinus1Contexts(
    const AlnGraph& graph,
    std::size_t start_node_id,
    std::size_t depth,
    std::unordered_set<std::size_t>& visited) const {

  std::vector<ContextInfo> contexts;
  const auto& successors = graph.outgoing(start_node_id);

  // If no successors, return empty context
  if (successors.empty()) {
    return {{"", 0}};
  }

  // For each immediate successor, collect context greedily
  // (take first available path, don't explore all branches for now to avoid explosion)
  for (std::size_t succ_id : successors) {
    // Cycle detection: skip if we've visited this node in current path
    if (visited.find(succ_id) != visited.end()) {
      continue;
    }

    // Greedily collect k-1 bases by following the first available path
    std::string collected_context;
    std::size_t current_node = succ_id;
    std::unordered_set<std::size_t> path_visited = visited;
    path_visited.insert(succ_id);

    // Collect up to 'depth' bases
    while (collected_context.size() < depth) {
      const AlnNode& curr_node = graph.node(current_node);

      // Extract bases from current node's sequence (up to what we need)
      std::size_t bases_needed = depth - collected_context.size();
      std::size_t bases_available = curr_node.sequence.size();
      std::size_t bases_to_take = std::min(bases_needed, bases_available);

      if (bases_to_take > 0) {
        collected_context += curr_node.sequence.substr(0, bases_to_take);
      }

      // If we've collected enough, break
      if (collected_context.size() >= depth) {
        break;
      }

      // Move to first available successor
      const auto& next_successors = graph.outgoing(current_node);
      if (next_successors.empty()) {
        break;  // No more successors, context incomplete
      }

      // Take first successor that hasn't been visited (greedy)
      bool found_next = false;
      for (std::size_t next_id : next_successors) {
        if (path_visited.find(next_id) == path_visited.end()) {
          current_node = next_id;
          path_visited.insert(next_id);
          found_next = true;
          break;
        }
      }

      if (!found_next) {
        break;  // All successors visited (cycle), stop
      }
    }

    // Add this context
    ContextInfo context_info;
    context_info.context = collected_context;
    context_info.successor_aln_id = succ_id;  // Track the immediate successor
    contexts.push_back(context_info);
  }

  return contexts;
}

// Expand uncovered nodes by creating F/R variants without context
void PathGuidedTransform::expandUncoveredNodes(
    AlnGraph& graph,
    const std::unordered_set<std::string>& uncovered_node_ids,
    NodeMapping& node_mapping,
    std::size_t k_minus_1,
    const piru::io::ImportedGraph& imported) {

  // Get imported node map
  std::unordered_map<std::string, const io::ImportedGraphNode*> imported_node_map;
  for (const auto& node : imported.nodes) {
    imported_node_map[node.id] = &node;
  }

  // For each uncovered node, create F/R variants without context
  for (const auto& uncovered_id : uncovered_node_ids) {
    const auto imported_it = imported_node_map.find(uncovered_id);
    if (imported_it == imported_node_map.end()) {
      continue;
    }
    const io::ImportedGraphNode& imported_node = *imported_it->second;

    // Create forward variant
    AlnNode fwd_node;
    fwd_node.label = uncovered_id + "_F";
    fwd_node.original_id = uncovered_id;
    fwd_node.is_reverse = false;
    fwd_node.sequence = imported_node.sequence;  // No context for uncovered nodes
    std::size_t fwd_id = graph.addNode(fwd_node);

    // Create reverse variant
    AlnNode rev_node;
    rev_node.label = uncovered_id + "_R";
    rev_node.original_id = uncovered_id;
    rev_node.is_reverse = true;
    rev_node.sequence = revcomp(imported_node.sequence);  // No context for uncovered nodes
    std::size_t rev_id = graph.addNode(rev_node);

    // Update node_mapping
    node_mapping[uncovered_id] = {fwd_id, rev_id};
  }
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
