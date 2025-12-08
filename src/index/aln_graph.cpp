// SPDX-License-Identifier: MIT
#include "index/aln_graph.hpp"

#include <algorithm>
#include <stdexcept>

namespace piru::index {

bool AlnGraph::validate() const {
    const std::size_t n = nodes_.size();
    bool ok = true;

    // Edge bounds and adjacency symmetry.
    for (std::size_t v = 0; v < n; ++v) {
        for (const auto u : out_edges_[v]) {
            if (u >= n) ok = false;
            const auto& incoming_u = in_edges_[u];
            const bool has_back = std::find(incoming_u.begin(), incoming_u.end(), v) != incoming_u.end();
            if (!has_back) ok = false;
        }
        for (const auto u : in_edges_[v]) {
            if (u >= n) ok = false;
            const auto& outgoing_u = out_edges_[u];
            const bool has_fwd = std::find(outgoing_u.begin(), outgoing_u.end(), v) != outgoing_u.end();
            if (!has_fwd) ok = false;
        }
    }

    // Optional metadata consistency: if chain_id present, linear_position should also be set.
    for (const auto& node : nodes_) {
        if (node.chain_id.has_value() != node.linear_position.has_value()) {
            ok = false;
        }
    }

    return ok;
}

}  // namespace piru::index
