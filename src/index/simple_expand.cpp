// SPDX-License-Identifier: MIT

#include "index/simple_expand.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "util/logging.hpp"

namespace piru::index {

namespace {

// Reverse complement a DNA sequence
std::string revcomp(const std::string& seq) {
  static const std::unordered_map<char, char> comp = {
      {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}, {'a', 't'},
      {'t', 'a'}, {'c', 'g'}, {'g', 'c'}, {'N', 'N'}, {'n', 'n'}};

  std::string rc;
  rc.reserve(seq.size());
  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    auto comp_it = comp.find(*it);
    rc += (comp_it != comp.end()) ? comp_it->second : 'N';
  }
  return rc;
}

}  // namespace

FlatGraph simpleExpandFlat(const piru::io::ImportedGraph& imported) {
  // Build mapping: original string ID -> index
  std::unordered_map<std::string, std::size_t> id_to_index;
  for (std::size_t i = 0; i < imported.nodes.size(); ++i) {
    id_to_index[imported.nodes[i].id] = i;
  }

  const std::size_t num_orig = imported.nodes.size();
  const std::uint32_t num_nodes = static_cast<std::uint32_t>(num_orig * 2);

  // --- Step 1: Build node arenas ---
  std::vector<char> seq_data;
  std::vector<char> name_data;
  std::vector<std::uint32_t> seq_offset(num_nodes);
  std::vector<std::uint32_t> seq_len(num_nodes);
  std::vector<std::uint32_t> name_offset(num_nodes);
  std::vector<std::uint16_t> name_len(num_nodes);
  std::vector<std::uint8_t> is_reverse(num_nodes);

  // Pre-estimate total sequence size
  std::size_t total_seq = 0;
  for (const auto& orig : imported.nodes) total_seq += orig.sequence.size() * 2;
  seq_data.reserve(total_seq);

  for (std::size_t i = 0; i < num_orig; ++i) {
    const auto& orig = imported.nodes[i];
    std::uint32_t fwd_id = static_cast<std::uint32_t>(forwardNodeId(i));
    std::uint32_t rev_id = static_cast<std::uint32_t>(reverseNodeId(i));

    // Forward node
    seq_offset[fwd_id] = static_cast<std::uint32_t>(seq_data.size());
    seq_len[fwd_id] = static_cast<std::uint32_t>(orig.sequence.size());
    seq_data.insert(seq_data.end(), orig.sequence.begin(), orig.sequence.end());

    name_offset[fwd_id] = static_cast<std::uint32_t>(name_data.size());
    name_len[fwd_id] = static_cast<std::uint16_t>(orig.id.size());
    name_data.insert(name_data.end(), orig.id.begin(), orig.id.end());
    is_reverse[fwd_id] = 0;

    // Reverse node
    std::string rc = revcomp(orig.sequence);
    seq_offset[rev_id] = static_cast<std::uint32_t>(seq_data.size());
    seq_len[rev_id] = static_cast<std::uint32_t>(rc.size());
    seq_data.insert(seq_data.end(), rc.begin(), rc.end());

    name_offset[rev_id] = static_cast<std::uint32_t>(name_data.size());
    name_len[rev_id] = static_cast<std::uint16_t>(orig.id.size());
    name_data.insert(name_data.end(), orig.id.begin(), orig.id.end());
    is_reverse[rev_id] = 1;
  }

  // --- Step 2: Build edges (CSR) ---
  // First pass: collect edges into adjacency lists
  std::vector<std::vector<std::uint32_t>> adj(num_nodes);

  for (const auto& orig_edge : imported.edges) {
    auto from_it = id_to_index.find(orig_edge.from);
    auto to_it = id_to_index.find(orig_edge.to);
    if (from_it == id_to_index.end() || to_it == id_to_index.end()) continue;

    std::size_t from_idx = from_it->second;
    std::size_t to_idx = to_it->second;

    // Forward edge
    std::uint32_t from_id = static_cast<std::uint32_t>(
        orig_edge.from_reverse ? reverseNodeId(from_idx) : forwardNodeId(from_idx));
    std::uint32_t to_id = static_cast<std::uint32_t>(
        orig_edge.to_reverse ? reverseNodeId(to_idx) : forwardNodeId(to_idx));
    adj[from_id].push_back(to_id);

    // Reverse complement edge
    std::uint32_t rev_from = static_cast<std::uint32_t>(
        orig_edge.to_reverse ? forwardNodeId(to_idx) : reverseNodeId(to_idx));
    std::uint32_t rev_to = static_cast<std::uint32_t>(
        orig_edge.from_reverse ? forwardNodeId(from_idx) : reverseNodeId(from_idx));
    adj[rev_from].push_back(rev_to);
  }

  // Dedup edges and build CSR
  std::vector<std::uint32_t> edge_target;
  std::vector<std::uint32_t> out_edge_offset(num_nodes + 1);
  for (std::uint32_t i = 0; i < num_nodes; ++i) {
    auto& edges = adj[i];
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    out_edge_offset[i] = static_cast<std::uint32_t>(edge_target.size());
    edge_target.insert(edge_target.end(), edges.begin(), edges.end());
  }
  out_edge_offset[num_nodes] = static_cast<std::uint32_t>(edge_target.size());
  adj.clear();

  // --- Step 3: Build paths ---
  const std::uint32_t num_paths = static_cast<std::uint32_t>(imported.paths.size() * 2);
  std::vector<std::uint32_t> step_data;
  std::vector<std::uint32_t> path_step_offset(num_paths + 1);
  std::vector<std::uint32_t> path_name_offset(num_paths);
  std::vector<std::uint16_t> path_name_len(num_paths);
  std::vector<std::uint64_t> path_length(num_paths, 0);

  for (std::size_t p = 0; p < imported.paths.size(); ++p) {
    const auto& orig_path = imported.paths[p];
    std::uint32_t fwd_p = static_cast<std::uint32_t>(p * 2);
    std::uint32_t rev_p = static_cast<std::uint32_t>(p * 2 + 1);

    // Forward path
    path_name_offset[fwd_p] = static_cast<std::uint32_t>(name_data.size());
    path_name_len[fwd_p] = static_cast<std::uint16_t>(orig_path.name.size());
    name_data.insert(name_data.end(), orig_path.name.begin(), orig_path.name.end());

    path_step_offset[fwd_p] = static_cast<std::uint32_t>(step_data.size());
    for (const auto& step : orig_path.steps) {
      auto it = id_to_index.find(step.segment_id);
      if (it == id_to_index.end()) continue;
      std::uint32_t nid = static_cast<std::uint32_t>(
          step.is_reverse ? reverseNodeId(it->second) : forwardNodeId(it->second));
      step_data.push_back(nid);
    }

    // Reverse path
    std::string rev_name = orig_path.name + "_reverse";
    path_name_offset[rev_p] = static_cast<std::uint32_t>(name_data.size());
    path_name_len[rev_p] = static_cast<std::uint16_t>(rev_name.size());
    name_data.insert(name_data.end(), rev_name.begin(), rev_name.end());

    path_step_offset[rev_p] = static_cast<std::uint32_t>(step_data.size());
    for (auto it = orig_path.steps.rbegin(); it != orig_path.steps.rend(); ++it) {
      auto idx_it = id_to_index.find(it->segment_id);
      if (idx_it == id_to_index.end()) continue;
      bool flipped = !it->is_reverse;
      std::uint32_t nid = static_cast<std::uint32_t>(
          flipped ? reverseNodeId(idx_it->second) : forwardNodeId(idx_it->second));
      step_data.push_back(nid);
    }
  }
  path_step_offset[num_paths] = static_cast<std::uint32_t>(step_data.size());

  return FlatGraph::fromRawArrays(
      num_nodes, num_paths,
      std::move(seq_data), std::move(seq_offset), std::move(seq_len),
      std::move(name_data), std::move(name_offset), std::move(name_len),
      std::move(is_reverse),
      std::move(edge_target), std::move(out_edge_offset),
      std::move(step_data), std::move(path_step_offset),
      std::move(path_name_offset), std::move(path_name_len), std::move(path_length));
}

}  // namespace piru::index
