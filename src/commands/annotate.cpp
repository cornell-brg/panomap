// SPDX-License-Identifier: MIT

#include "commands/annotate.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/parse.hpp"
#include "index/aln_graph.hpp"
#include "index/linearizer.hpp"
#include "index/simple_expand.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/regions/bed_parser.hpp"
#include "util/logging.hpp"

namespace {

/// Result of walking a path interval to find overlapping nodes.
struct IntervalWalk {
    std::vector<std::size_t> node_ids;  // AlnGraph node IDs in path order
    std::size_t start_offset{0};        // offset into first node where BED starts
    std::size_t end_offset{0};          // offset into last node where BED ends
    std::int64_t interval_length{0};    // bed_end - bed_start
};

/// Walk a path's steps to find nodes overlapping [bed_start, bed_end).
/// Returns the node IDs, start/end offsets into boundary nodes, and interval length.
IntervalWalk walkInterval(const piru::index::AlnGraph& graph,
                          const piru::index::AlnPath& path,
                          std::int64_t bed_start, std::int64_t bed_end) {
    IntervalWalk result;
    result.interval_length = bed_end - bed_start;

    std::int64_t offset = 0;
    for (const auto& step : path.steps) {
        std::size_t node_id = std::stoull(step.node_id);
        if (node_id >= graph.nodeCount()) continue;

        std::int64_t node_len = static_cast<std::int64_t>(graph.node(node_id).sequence.size());
        std::int64_t node_start = offset;
        std::int64_t node_end = offset + node_len;

        // Check overlap: [node_start, node_end) ^ [bed_start, bed_end)
        if (node_end > bed_start && node_start < bed_end) {
            if (result.node_ids.empty()) {
                result.start_offset = static_cast<std::size_t>(bed_start - node_start);
            }
            result.node_ids.push_back(node_id);
            // Update end_offset for each node -last one wins
            std::int64_t this_node_bed_end = std::min(bed_end, node_end);
            result.end_offset = static_cast<std::size_t>(this_node_bed_end - node_start);
        }

        // Past the interval -stop early
        if (node_start >= bed_end) break;

        offset = node_end;
    }

    return result;
}


/// Result of the mini DP: the equivalent interval on a haplotype path.
struct EquivalentInterval {
    std::vector<std::size_t> node_ids;  // selected node IDs in coord order
    std::int64_t start_coord{0};        // coord of first selected node
    std::int64_t end_coord{0};          // coord of last selected node + its length
    std::size_t chain_length{0};        // number of nodes in the chain
};

/// Mini DP: find the longest increasing subsequence of (node_id, coord) candidates
/// on a haplotype path, with a max-gap constraint between consecutive entries.
/// Candidates are (node_id, coord) pairs from step 2's lookup.
EquivalentInterval findEquivalentInterval(
        const std::vector<std::pair<std::size_t, std::int64_t>>& candidates,
        const piru::index::AlnGraph& graph,
        std::int64_t max_gap) {
    if (candidates.empty()) return {};

    // Sort by coord
    auto sorted = candidates;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    const std::size_t n = sorted.size();
    std::vector<std::size_t> dp_len(n, 1);        // longest chain ending at i
    std::vector<std::int64_t> dp_start(n);         // start coord of chain ending at i
    std::vector<std::size_t> pred(n, n);           // predecessor index (n = no pred)

    // Base case: each candidate starts its own chain
    for (std::size_t i = 0; i < n; ++i) {
        dp_start[i] = sorted[i].second;
    }

    // DP: maximize chain length, tiebreak by minimum span
    for (std::size_t i = 1; i < n; ++i) {
        for (std::size_t j = i; j-- > 0;) {
            std::int64_t gap = sorted[i].second - sorted[j].second;
            if (gap > max_gap) break;  // sorted by coord, so all earlier j are further away

            std::size_t new_len = dp_len[j] + 1;
            if (new_len > dp_len[i]) {
                dp_len[i] = new_len;
                dp_start[i] = dp_start[j];
                pred[i] = j;
            } else if (new_len == dp_len[i] && dp_start[j] > dp_start[i]) {
                // Same length but tighter span (later start = smaller span)
                dp_start[i] = dp_start[j];
                pred[i] = j;
            }
        }
    }

    // Find best endpoint: max length, then min span
    std::size_t best_idx = 0;
    for (std::size_t i = 1; i < n; ++i) {
        if (dp_len[i] > dp_len[best_idx]) {
            best_idx = i;
        } else if (dp_len[i] == dp_len[best_idx]) {
            std::int64_t span_i = sorted[i].second - dp_start[i];
            std::int64_t span_best = sorted[best_idx].second - dp_start[best_idx];
            if (span_i < span_best) best_idx = i;
        }
    }

    // Backtrack to get chain
    EquivalentInterval result;
    result.chain_length = dp_len[best_idx];

    std::vector<std::size_t> chain_indices;
    for (std::size_t idx = best_idx; idx < n; idx = pred[idx]) {
        chain_indices.push_back(idx);
    }
    std::reverse(chain_indices.begin(), chain_indices.end());

    for (std::size_t idx : chain_indices) {
        result.node_ids.push_back(sorted[idx].first);
    }

    // Interval bounds
    std::size_t first = chain_indices.front();
    std::size_t last = chain_indices.back();
    result.start_coord = sorted[first].second;
    result.end_coord = sorted[last].second +
        static_cast<std::int64_t>(graph.node(sorted[last].first).sequence.size());

    return result;
}

}  // namespace

int handle_annotate(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru annotate [options] <graph-file>";
    config.positional_help = {"<graph-file>      GFA graph file"};
    config.options = {
        {'h', "help", false, "Show this help message"},
        {'b', "bed", true, "BED file with target regions (path_name start end)"},
        {'o', "output", true, "Output annotation file (.pira)"},
        {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
        {0, "original-ids", false, "Output original GFA node IDs instead of AlnGraph IDs"},
    };
    config.on_error = [](const std::string&) {
        std::cerr << "annotate: invalid option\n";
    };

    if (!piru::cli::parse_args(args, config, parsed)) {
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.values.count("help")) {
        piru::cli::print_help(config, std::cout);
        return 0;
    }

    const bool verbose = parsed.values.count("verbose") > 0;
    if (verbose) {
        piru::logger.set_level(piru::LogLevel::DEBUG);
    }
    const bool use_original_ids = parsed.values.count("original-ids") > 0;

    // --- Validate required args ---
    if (parsed.positionals.empty()) {
        LOG_ERROR("annotate: missing required <graph-file>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    const std::string graph_path = parsed.positionals.front();

    if (!parsed.values.count("bed")) {
        LOG_ERROR("annotate: missing required --bed <file>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    const std::string bed_path = parsed.values.at("bed");

    const std::string output_path =
        parsed.values.count("output") ? parsed.values.at("output") : "annotate.pira";

    // --- Parse BED ---
    const auto bed_records = piru::io::parse_bed(bed_path);
    if (bed_records.empty()) {
        LOG_ERROR("annotate: no valid BED records found");
        return 1;
    }

    // --- Load GFA ---
    auto loader = piru::io::make_graph_loader(graph_path);
    if (!loader) {
        LOG_ERROR("annotate: unsupported graph format for '" + graph_path + "'");
        return 1;
    }

    piru::io::ImportedGraph imported;
    imported.flavor = piru::io::ImportedGraphFlavor::kVg;
    if (!loader->load(imported)) {
        LOG_ERROR("annotate: failed to read graph file '" + graph_path + "'");
        return 1;
    }
    LOG_INFO("Loaded graph: " + std::to_string(imported.nodes.size()) + " nodes, " +
             std::to_string(imported.paths.size()) + " paths");

    // --- +/-expand -> AlnGraph ---
    auto aln_graph = piru::index::simpleExpand(imported);
    LOG_INFO("AlnGraph: " + std::to_string(aln_graph.nodeCount()) + " nodes (+/-expand), " +
             std::to_string(aln_graph.pathCount()) + " paths");

    // --- Build linearization coords (path walks + reverse map) ---
    std::vector<std::vector<piru::index::LinearCoordinate>> lin_coords(aln_graph.nodeCount());
    std::unordered_map<std::string, std::size_t> path_name_to_id;

    for (std::size_t path_idx = 0; path_idx < aln_graph.pathCount(); ++path_idx) {
        const auto& path = aln_graph.paths()[path_idx];
        path_name_to_id[path.name] = path_idx;

        std::int64_t offset = 0;
        for (const auto& step : path.steps) {
            std::size_t node_id = std::stoull(step.node_id);
            if (node_id >= aln_graph.nodeCount()) continue;

            lin_coords[node_id].emplace_back(path_idx, offset);
            offset += static_cast<std::int64_t>(aln_graph.node(node_id).sequence.size());
        }
    }

    // --- Step 1: Project BED intervals onto reference path ---
    for (const auto& rec : bed_records) {
        auto it = path_name_to_id.find(rec.path_name);
        if (it == path_name_to_id.end()) {
            LOG_WARN("BED path '" + rec.path_name + "' not found in graph, skipping");
            continue;
        }
        std::size_t path_idx = it->second;
        const auto& path = aln_graph.paths()[path_idx];

        auto walk = walkInterval(aln_graph, path, rec.start, rec.end);

        LOG_INFO("Step 1: " + rec.path_name + ":" + std::to_string(rec.start) + "-" +
                 std::to_string(rec.end) + " -> " + std::to_string(walk.node_ids.size()) +
                 " nodes, start_offset=" + std::to_string(walk.start_offset) +
                 ", end_offset=" + std::to_string(walk.end_offset));

        // --- Step 2: For each ROI node, look up coords on all other forward paths ---
        // Group by path: path_id -> vector of (node_id, coord)
        std::map<std::size_t, std::vector<std::pair<std::size_t, std::int64_t>>> path_candidates;

        for (std::size_t node_id : walk.node_ids) {
            for (const auto& lc : lin_coords[node_id]) {
                // Skip the reference path itself
                if (lc.path_id == path_idx) continue;
                // Skip reverse paths (odd indices) -BED is forward-strand
                if (lc.path_id % 2 != 0) continue;

                path_candidates[lc.path_id].emplace_back(node_id, lc.ref_coord);
            }
        }

        LOG_INFO("Step 2: ROI nodes appear on " + std::to_string(path_candidates.size()) +
                 " other forward paths");

        // --- Step 3: Mini DP per haplotype path to find equivalent interval ---
        std::int64_t max_gap = walk.interval_length;  // no consecutive pair further than interval length

        LOG_INFO("Step 3: ref interval " + rec.path_name + ":[" +
                 std::to_string(rec.start) + ", " + std::to_string(rec.end) +
                 ") L=" + std::to_string(walk.interval_length) +
                 ", " + std::to_string(walk.node_ids.size()) + " ref nodes");

        // --- Steps 4+5: Walk equivalent intervals, union all node sets ---
        std::set<std::size_t> roi_nodes;

        // Add reference walk nodes
        for (std::size_t nid : walk.node_ids) {
            roi_nodes.insert(nid);
        }
        std::size_t ref_count = roi_nodes.size();

        for (const auto& [pid, candidates] : path_candidates) {
            auto equiv = findEquivalentInterval(candidates, aln_graph, max_gap);

            std::int64_t span = equiv.end_coord - equiv.start_coord;
            LOG_INFO("  " + aln_graph.paths()[pid].name +
                     ": [" + std::to_string(equiv.start_coord) + ", " +
                     std::to_string(equiv.end_coord) + ") span=" + std::to_string(span) +
                     ", " + std::to_string(equiv.chain_length) + "/" +
                     std::to_string(candidates.size()) + " candidates kept");

            // Step 4: Walk this path's equivalent interval to get full node list
            const auto& hap_path = aln_graph.paths()[pid];
            auto hap_walk = walkInterval(aln_graph, hap_path,
                                         equiv.start_coord, equiv.end_coord);

            // Step 5: Add to union
            for (std::size_t nid : hap_walk.node_ids) {
                roi_nodes.insert(nid);
            }
        }

        // Add reverse counterparts for all ROI nodes (+/-expand)
        std::set<std::size_t> roi_nodes_with_rev;
        for (std::size_t nid : roi_nodes) {
            roi_nodes_with_rev.insert(nid);
            std::size_t orig = piru::index::originalIndex(nid);
            if (piru::index::isReverseNode(nid)) {
                roi_nodes_with_rev.insert(piru::index::forwardNodeId(orig));
            } else {
                roi_nodes_with_rev.insert(piru::index::reverseNodeId(orig));
            }
        }

        LOG_INFO("Step 4+5: ROI union = " + std::to_string(roi_nodes.size()) +
                 " fwd nodes (" + std::to_string(ref_count) + " from ref, " +
                 std::to_string(roi_nodes.size() - ref_count) + " from haplotypes), " +
                 std::to_string(roi_nodes_with_rev.size()) + " total with +/-rev");

        // --- Output annotation ---
        std::ofstream ofs(output_path, std::ios::app);
        if (!ofs.is_open()) {
            LOG_ERROR("annotate: cannot open output file '" + output_path + "'");
            return 1;
        }

        ofs << rec.path_name << ":" << rec.start << "-" << rec.end;
        if (use_original_ids) {
            // Deduplicate to original GFA node names (forward/reverse collapse to same original)
            std::set<std::string> orig_names;
            for (std::size_t nid : roi_nodes_with_rev) {
                std::size_t idx = piru::index::originalIndex(nid);
                orig_names.insert(imported.nodes[idx].id);
            }
            bool first = true;
            for (const auto& name : orig_names) {
                ofs << (first ? " " : ",") << name;
                first = false;
            }
        } else {
            bool first = true;
            for (std::size_t nid : roi_nodes_with_rev) {
                ofs << (first ? " " : ",") << nid;
                first = false;
            }
        }
        ofs << "\n";
    }

    LOG_INFO("Annotation written to " + output_path);

    return 0;
}
