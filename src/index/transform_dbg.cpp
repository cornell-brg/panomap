// SPDX-License-Identifier: MIT
#include "index/transform_dbg.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "util/logging.hpp"

namespace piru::index {
namespace {

std::string revcomp(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        const char c = *it;
        switch (c) {
            case 'A':
                rc.push_back('T');
                break;
            case 'C':
                rc.push_back('G');
                break;
            case 'G':
                rc.push_back('C');
                break;
            case 'T':
                rc.push_back('A');
                break;
            case 'a':
                rc.push_back('t');
                break;
            case 'c':
                rc.push_back('g');
                break;
            case 'g':
                rc.push_back('c');
                break;
            case 't':
                rc.push_back('a');
                break;
            default:
                rc.push_back('N');
                break;
        }
    }
    return rc;
}

}  // namespace

AlnGraph transformDbg(const io::ImportedGraph& imported, std::size_t graph_k, std::size_t pore_k) {
    // Validation: pore model k cannot exceed DBG k
    if (pore_k > graph_k) {
        throw std::invalid_argument(
            "pore_k (" + std::to_string(pore_k) +
            ") cannot be greater than graph_k (" + std::to_string(graph_k) + ")");
    }

    // Calculate how much to trim from node sequences
    // DBG nodes have (graph_k - 1) overlap, but we only need (pore_k - 1)
    // So trim the last k_delta = (graph_k - pore_k) bases
    const std::size_t k_delta = graph_k - pore_k;

    AlnGraph graph;
    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> node_lookup;
    node_lookup.reserve(imported.nodes.size());

    // Create directional nodes (forward and reverse) with trimmed sequences
    for (const auto& n : imported.nodes) {
        // Check node is long enough after trimming
        if (n.sequence.size() <= k_delta) {
            // Node too short - skip with warning
            std::ostringstream msg;
            msg << "Skipping node '" << n.id
                << "' (length " << n.sequence.size()
                << ") - too short after trimming k_delta=" << k_delta
                << " (graph_k=" << graph_k << ", pore_k=" << pore_k << ")";
            LOG_WARN(msg.str());
            continue;
        }

        // Trim last k_delta bases to preserve only pore_k-1 overlap
        const std::string trimmed_seq = n.sequence.substr(0, n.sequence.size() - k_delta);

        AlnNode fwd;
        fwd.label = n.id;
        fwd.original_id = n.id;
        fwd.is_reverse = false;
        fwd.sequence = trimmed_seq;
        const std::size_t fwd_id = graph.addNode(std::move(fwd));

        AlnNode rev;
        rev.label = n.id;
        rev.original_id = n.id;
        rev.is_reverse = true;
        rev.sequence = revcomp(trimmed_seq);  // Reverse complement of trimmed sequence
        const std::size_t rev_id = graph.addNode(std::move(rev));

        node_lookup[n.id] = {fwd_id, rev_id};
    }

    // Add edges respecting orientations with adjusted overlaps
    for (const auto& e : imported.edges) {
        const auto it_from = node_lookup.find(e.from);
        const auto it_to = node_lookup.find(e.to);
        if (it_from == node_lookup.end() || it_to == node_lookup.end()) {
            continue;
        }

        const auto from_id = e.from_reverse ? it_from->second.second : it_from->second.first;
        const auto to_id = e.to_reverse ? it_to->second.second : it_to->second.first;

        // Adjust overlap: original DBG overlap is (graph_k - 1), after trimming it's (pore_k - 1)
        std::size_t overlap = e.overlap_bases.value_or(graph_k - 1);
        if (overlap >= k_delta) {
            overlap -= k_delta;
        } else {
            // Safety: if overlap is somehow less than k_delta, set to 0
            overlap = 0;
        }

        graph.addEdge({from_id, to_id, overlap});
    }

    return graph;
}

}  // namespace piru::index
