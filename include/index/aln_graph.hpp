// SPDX-License-Identifier: MIT
// Directional alignment graph representation used during indexing.

#pragma once

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace piru::index {

struct AlnPathStep {
    std::string node_id;
    bool is_reverse{false}; // Indicates if this step traverses the node in reverse
};

struct AlnPath {
    std::string name;
    std::vector<AlnPathStep> steps;
    // Overlap length in bases for each step, if explicitly tracked (e.g., from ImportedPath).
    // The size would be steps.size() - 1 for overlaps between steps.
    std::vector<std::size_t> overlaps;
};

// ... (rest of the file)

struct AlnNode {
    std::size_t id{0};                // Sequential ID within the graph.
    std::string label;                // Friendly label (default: original node id).
    std::string original_id;          // Original node id from the imported graph.
    bool is_reverse{false};           // True if this node is the reverse complement.
    std::string sequence;             // Sequence including k-1 overlaps where applicable.
    std::optional<std::int64_t> chain_id;
    std::optional<std::int64_t> linear_position;
};

struct AlnEdge {
    std::size_t from{0};
    std::size_t to{0};
    std::size_t overlap_bases{0};
};

class AlnGraph {
public:
    std::size_t addNode(AlnNode node) {
        const std::size_t idx = nodes_.size();
        node.id = idx;
        nodes_.push_back(std::move(node));
        out_edges_.emplace_back();
        in_edges_.emplace_back();
        return idx;
    }

    void addEdge(AlnEdge edge) {
        if (edge.from >= nodes_.size() || edge.to >= nodes_.size()) return;
        edges_.push_back(edge);
        out_edges_[edge.from].push_back(edge.to);
        in_edges_[edge.to].push_back(edge.from);
    }

    void clear() {
        nodes_.clear();
        edges_.clear();
        out_edges_.clear();
        in_edges_.clear();
    }

    std::size_t nodeCount() const { return nodes_.size(); }
    std::size_t edgeCount() const { return edges_.size(); }

    const AlnNode& node(std::size_t idx) const { return nodes_.at(idx); }
    AlnNode& mutableNode(std::size_t idx) { return nodes_.at(idx); }

    const std::vector<std::size_t>& outgoing(std::size_t idx) const { return out_edges_.at(idx); }
    const std::vector<std::size_t>& incoming(std::size_t idx) const { return in_edges_.at(idx); }

    // O(degree) edge existence check using adjacency list
    bool hasEdge(std::size_t from, std::size_t to) const {
        if (from >= out_edges_.size()) return false;
        const auto& neighbors = out_edges_[from];
        return std::find(neighbors.begin(), neighbors.end(), to) != neighbors.end();
    }

    const std::vector<AlnEdge>& edges() const { return edges_; }

    void addPath(AlnPath path) { paths_.push_back(std::move(path)); } // Fixed std::std::move
    const std::vector<AlnPath>& paths() const { return paths_; }
    std::size_t pathCount() const { return paths_.size(); }

    // Basic consistency checks: edge bounds, adjacency symmetry, optional metadata presence.
    bool validate() const;

private:
    std::vector<AlnNode> nodes_;
    std::vector<AlnEdge> edges_;
    std::vector<std::vector<std::size_t>> out_edges_;
    std::vector<std::vector<std::size_t>> in_edges_;
    std::vector<AlnPath> paths_;
};

} // namespace piru::index
