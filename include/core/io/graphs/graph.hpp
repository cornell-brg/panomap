// SPDX-License-Identifier: MIT
// Minimal in-memory representation of graphs parsed from files.

#pragma once

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace piru::io {

struct ImportedGraphNode {
  std::string id;
  std::string sequence;
};

struct ImportedGraphEdge {
  std::string from;
  std::string to;
  bool from_reverse{false};
  bool to_reverse{false};
  // Raw overlap: CIGAR for GFA links or stringified length for vg.
  // Defaults to "0" when unspecified.
  std::string overlap{"0"};
  // Parsed overlap length when available (GFA CIGAR parsed or vg overlap length).
  std::optional<std::size_t> overlap_bases;
};

struct ImportedPathStep {
  std::string segment_id;
  bool is_reverse{false};
};

struct ImportedPath {
  std::string name;
  std::vector<ImportedPathStep> steps;
  // Optional per-edge overlaps as raw CIGAR strings; size is steps.size() - 1 when present.
  std::vector<std::string> overlaps;
};

class ImportedGraph {
public:
  void add_node(ImportedGraphNode node) {
    lookup_[node.id] = nodes.size();
    nodes.push_back(std::move(node));
  }

  void add_edge(ImportedGraphEdge edge) { edges.push_back(std::move(edge)); }

  void add_path(ImportedPath path) { paths.push_back(std::move(path)); }

  const ImportedGraphNode* find_node(const std::string& id) const {
    const auto it = lookup_.find(id);
    if (it == lookup_.end()) return nullptr;
    return &nodes[it->second];
  }

  void clear() {
    nodes.clear();
    edges.clear();
    lookup_.clear();
  }

  std::vector<ImportedGraphNode> nodes;
  std::vector<ImportedGraphEdge> edges;
  std::vector<ImportedPath> paths;

private:
  std::unordered_map<std::string, std::size_t> lookup_;
};

}  // namespace piru::io
