// SPDX-License-Identifier: MIT
#include "index/aln_graph.hpp"

#include <algorithm>
#include <unordered_set>

namespace piru::index {

bool AlnGraph::validate() const {
    const std::size_t n = nodes_.size();
    bool ok = true;

    // Edge bounds and adjacency symmetry.
    for (std::size_t v = 0; v < n; ++v) {
        for (const auto u : out_edges_[v]) {
            if (u >= n) {
                ok = false;
                break;  // Exit early if an invalid edge is found
            }
            const auto& incoming_u = in_edges_[u];
            const bool has_back =
                std::find(incoming_u.begin(), incoming_u.end(), v) != incoming_u.end();
            if (!has_back) {
                ok = false;
                break;  // Exit early if asymmetry is found
            }
        }
        if (!ok) break;  // Exit outer loop if an issue was found
        for (const auto u : in_edges_[v]) {
            if (u >= n) {
                ok = false;
                break;  // Exit early if an invalid edge is found
            }
            const auto& outgoing_u = out_edges_[u];
            const bool has_fwd =
                std::find(outgoing_u.begin(), outgoing_u.end(), v) != outgoing_u.end();
            if (!has_fwd) {
                ok = false;
                break;  // Exit early if asymmetry is found
            }
        }
        if (!ok) break;  // Exit outer loop if an issue was found
    }

    if (!ok) return false;  // If basic graph structure is invalid, no need to check further

    // Optional metadata consistency: if chain_id present, linear_position should also be set.
    for (const auto& node : nodes_) {
        if (node.chain_id.has_value() != node.linear_position.has_value()) {
            ok = false;
            break;
        }
    }

    if (!ok) return false;

    // Build label set for O(1) lookups (instead of O(N) per step)
    std::unordered_set<std::string> label_set;
    label_set.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        label_set.insert(node.label);
    }

    // Path node references and overlap semantics.
    for (const auto& path : paths_) {
        // Validate overlap vector size
        if (!path.overlaps.empty() && path.overlaps.size() != path.steps.size() - 1) {
            ok = false;
            break;
        }

        for (const auto& step : path.steps) {
            // O(1) lookup instead of O(N) scan
            if (label_set.find(step.node_id) == label_set.end()) {
                ok = false;
                break;
            }
        }
        if (!ok) break;
    }

    return ok;
}

}  // namespace piru::index
