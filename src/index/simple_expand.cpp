// SPDX-License-Identifier: MIT

#include "index/simple_expand.hpp"

#include <algorithm>
#include <unordered_map>

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

AlnGraph simpleExpand(const piru::io::ImportedGraph& imported) {
    AlnGraph graph;

    // Build mapping: original string ID → index (0-based)
    // Needed because edges/paths reference nodes by string ID, but we need the
    // numeric index to compute AlnGraph node IDs (index*2 for fwd, index*2+1 for rev)
    std::unordered_map<std::string, std::size_t> id_to_index;
    for (std::size_t i = 0; i < imported.nodes.size(); ++i) {
        id_to_index[imported.nodes[i].id] = i;
    }

    // -------------------------------------------------------------------------
    // Step 1: Split each bidirectional node to forward + reverse nodes
    // -------------------------------------------------------------------------
    // Node ID scheme:
    // - fwd: original_index * 2
    // - rev: original_index * 2 + 1

    const std::size_t num_nodes = imported.nodes.size() * 2;
    graph.reserveNodes(num_nodes);

    for (std::size_t i = 0; i < imported.nodes.size(); ++i) {
        const auto& orig = imported.nodes[i];

        // Forward node at ID = i * 2
        AlnNode fwd;
        fwd.label = orig.id + "+";  // TODO: remove later, not needed
        fwd.original_id = orig.id;  // TODO: remove later, not needed
        fwd.is_reverse = false;
        fwd.sequence = orig.sequence;
        graph.setNode(forwardNodeId(i), fwd);

        // Reverse node at ID = i * 2 + 1
        AlnNode rev;
        rev.label = orig.id + "-";  // TODO: remove later, not needed
        rev.original_id = orig.id;  // TODO: remove later, not needed
        rev.is_reverse = true;
        rev.sequence = revcomp(orig.sequence);
        graph.setNode(reverseNodeId(i), rev);
    }

    // -------------------------------------------------------------------------
    // Step 2: Transform edges (handle orientation flags)
    // -------------------------------------------------------------------------
    // ImportedGraphEdge has from_reverse and to_reverse flags
    // We need to create edges between the appropriate ±nodes

    for (const auto& orig_edge : imported.edges) {
        auto from_it = id_to_index.find(orig_edge.from);
        auto to_it = id_to_index.find(orig_edge.to);
        if (from_it == id_to_index.end() || to_it == id_to_index.end()) {
            continue;  // Skip edges with unknown nodes
        }

        std::size_t from_idx = from_it->second;
        std::size_t to_idx = to_it->second;

        // Determine which ±node to use based on orientation
        std::size_t from_aln_id =
            orig_edge.from_reverse ? reverseNodeId(from_idx) : forwardNodeId(from_idx);
        std::size_t to_aln_id =
            orig_edge.to_reverse ? reverseNodeId(to_idx) : forwardNodeId(to_idx);

        // Add forward edge
        if (!graph.hasEdge(from_aln_id, to_aln_id)) {
            AlnEdge edge;
            edge.from = from_aln_id;
            edge.to = to_aln_id;
            edge.overlap_bases = orig_edge.overlap_bases.value_or(0);
            graph.addEdge(edge);
        }

        // Add reverse complement edge: rev(to) → rev(from)
        // If original edge is A+ → B+, reverse is B- → A-
        std::size_t rev_from_aln_id =
            orig_edge.to_reverse ? forwardNodeId(to_idx) : reverseNodeId(to_idx);
        std::size_t rev_to_aln_id =
            orig_edge.from_reverse ? forwardNodeId(from_idx) : reverseNodeId(from_idx);

        if (!graph.hasEdge(rev_from_aln_id, rev_to_aln_id)) {
            AlnEdge rev_edge;
            rev_edge.from = rev_from_aln_id;
            rev_edge.to = rev_to_aln_id;
            rev_edge.overlap_bases = orig_edge.overlap_bases.value_or(0);
            graph.addEdge(rev_edge);
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: Transform paths (create forward + reverse paths)
    // -------------------------------------------------------------------------

    for (const auto& orig_path : imported.paths) {
        // Forward path
        AlnPath fwd_path;
        fwd_path.name = orig_path.name;
        for (const auto& step : orig_path.steps) {
            auto it = id_to_index.find(step.segment_id);
            if (it == id_to_index.end()) continue;

            std::size_t aln_id =
                step.is_reverse ? reverseNodeId(it->second) : forwardNodeId(it->second);
            AlnPathStep aln_step;
            aln_step.node_id = std::to_string(aln_id);  // Use numeric AlnGraph node ID
            aln_step.is_reverse = false;  // Always + since orientation is encoded in node choice
            fwd_path.steps.push_back(aln_step);
        }
        graph.addPath(fwd_path);

        // Reverse path (walk in reverse order with flipped orientations)
        AlnPath rev_path;
        rev_path.name = orig_path.name + "_reverse";
        for (auto it = orig_path.steps.rbegin(); it != orig_path.steps.rend(); ++it) {
            const auto& step = *it;
            auto idx_it = id_to_index.find(step.segment_id);
            if (idx_it == id_to_index.end()) continue;

            // Flip orientation for reverse path
            bool flipped = !step.is_reverse;
            std::size_t aln_id =
                flipped ? reverseNodeId(idx_it->second) : forwardNodeId(idx_it->second);
            AlnPathStep aln_step;
            aln_step.node_id = std::to_string(aln_id);  // Use numeric AlnGraph node ID
            aln_step.is_reverse = false;  // Always + since orientation is encoded in node choice
            rev_path.steps.push_back(aln_step);
        }
        graph.addPath(rev_path);
    }

    LOG_DEBUG("ImportedGraph to AlnGraph transformation done.");
    return graph;
}

}  // namespace piru::index
