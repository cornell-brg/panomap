/**
 * annotate.cpp
 *
 * Projects BED intervals onto graph node space and 1D canonical coordinates.
 * Two modes:
 *   strict: exact nodes from the BED interval's haplotype path
 *   union:  all nodes whose 1D coord falls in the interval range
 *
 * Output: .pira v2 format with 1D intervals.
 *
 * Related:
 *  - annotate.hpp
 *  - sort_1d.hpp (1D coordinate import)
 *  - pira_parser.hpp (reads .pira files)
 *
 * SPDX-License-Identifier: MIT
 */

#include "commands/annotate.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/parse.hpp"
#include "index/aln_graph.hpp"
#include "index/linearizer.hpp"
#include "index/simple_expand.hpp"
#include "index/sort_1d.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/regions/bed_parser.hpp"
#include "util/logging.hpp"

namespace {

// Walk a path to collect node IDs overlapping [bed_start, bed_end).
std::vector<std::size_t> walkInterval(const piru::index::AlnGraph& graph,
                                       const piru::index::AlnPath& path,
                                       std::int64_t bed_start, std::int64_t bed_end) {
  std::vector<std::size_t> node_ids;
  std::int64_t offset = 0;
  for (const auto& step : path.steps) {
    std::size_t node_id = std::stoull(step.node_id);
    if (node_id >= graph.nodeCount()) continue;

    std::int64_t node_len = static_cast<std::int64_t>(graph.node(node_id).sequence.size());
    std::int64_t node_start = offset;
    std::int64_t node_end = offset + node_len;

    if (node_end > bed_start && node_start < bed_end) {
      node_ids.push_back(node_id);
    }
    if (node_start >= bed_end) break;
    offset = node_end;
  }
  return node_ids;
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
      {'o', "output", true, "Output annotation file (.pira, default: annotate.pira)"},
      {'\0', "1d-coords-file", true, "1D coords TSV (from odgi sort --path-sgd-layout)"},
      {'\0', "with-nodes", false, "Include node set in output (for --roi usage)"},
      {'\0', "mode", true, "Node selection mode with --with-nodes: strict (default) or union"},
      {'v', "verbose", false, "Enable verbose logging"},
  };
  config.on_error = [](const std::string&) { std::cerr << "annotate: invalid option\n"; };

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
  const bool with_nodes = parsed.values.count("with-nodes") > 0;
  const std::string mode =
      parsed.values.count("mode") ? parsed.values.at("mode") : "strict";

  if (mode != "strict" && mode != "union") {
    LOG_ERROR("annotate: --mode must be 'strict' or 'union'");
    return 1;
  }

  /* Parse BED */
  const auto bed_records = piru::io::parse_bed(bed_path);
  if (bed_records.empty()) {
    LOG_ERROR("annotate: no valid BED records found");
    return 1;
  }

  /* Load GFA */
  auto loader = piru::io::make_graph_loader(graph_path);
  if (!loader) {
    LOG_ERROR("annotate: unsupported graph format for '" + graph_path + "'");
    return 1;
  }
  piru::io::ImportedGraph imported;
  if (!loader->load(imported)) {
    LOG_ERROR("annotate: failed to read graph file '" + graph_path + "'");
    return 1;
  }
  LOG_INFO("Loaded graph: " + std::to_string(imported.nodes.size()) + " nodes, " +
           std::to_string(imported.paths.size()) + " paths");

  auto aln_graph = piru::index::simpleExpand(imported);
  LOG_INFO("AlnGraph: " + std::to_string(aln_graph.nodeCount()) + " nodes, " +
           std::to_string(aln_graph.pathCount()) + " paths");

  /* Build path name -> id map */
  std::unordered_map<std::string, std::size_t> path_name_to_id;
  for (std::size_t i = 0; i < aln_graph.pathCount(); ++i) {
    path_name_to_id[aln_graph.paths()[i].name] = i;
  }

  /* Load 1D coords if provided */
  std::vector<float> node_1d_coords;
  bool has_1d = false;
  if (parsed.values.count("1d-coords-file")) {
    node_1d_coords = piru::index::import_1d_coords_odgi(
        parsed.values.at("1d-coords-file"), aln_graph.nodeCount());
    has_1d = true;
  }

  /* Write .pira v2 */
  std::ofstream ofs(output_path);
  if (!ofs.is_open()) {
    LOG_ERROR("annotate: cannot open output file '" + output_path + "'");
    return 1;
  }
  // TODO: set to v1.0 at release
  ofs << "#pira v2" << (has_1d ? " coords=1d" : "") << (with_nodes ? " nodes=" + mode : "") << "\n";

  for (const auto& rec : bed_records) {
    auto it = path_name_to_id.find(rec.path_name);
    if (it == path_name_to_id.end()) {
      // Try with + suffix (simpleExpand appends +/- to path names)
      it = path_name_to_id.find(rec.path_name + "+");
      if (it == path_name_to_id.end()) {
        LOG_WARN("BED path '" + rec.path_name + "' not found in graph, skipping");
        continue;
      }
    }
    std::size_t path_idx = it->second;
    const auto& path = aln_graph.paths()[path_idx];

    /* Walk the path to get nodes in BED interval */
    auto walk_nodes = walkInterval(aln_graph, path, rec.start, rec.end);
    if (walk_nodes.empty()) {
      LOG_WARN("No nodes found for " + rec.path_name + ":" +
               std::to_string(rec.start) + "-" + std::to_string(rec.end));
      continue;
    }

    /* Compute 1D interval from walked nodes */
    double coord_1d_start = 0.0, coord_1d_end = 0.0;
    if (has_1d) {
      coord_1d_start = std::numeric_limits<double>::max();
      coord_1d_end = std::numeric_limits<double>::lowest();
      for (std::size_t nid : walk_nodes) {
        if (nid < node_1d_coords.size()) {
          double node_start = node_1d_coords[nid];
          double node_end = node_start + aln_graph.node(nid).sequence.size();
          coord_1d_start = std::min(coord_1d_start, node_start);
          coord_1d_end = std::max(coord_1d_end, node_end);
        }
      }
    }

    /* Collect output nodes (only if --with-nodes) */
    std::set<std::size_t> roi_with_rev;

    if (with_nodes) {
      std::set<std::size_t> roi_nodes;

      if (mode == "strict") {
        for (std::size_t nid : walk_nodes) {
          roi_nodes.insert(nid);
        }
      } else {
        if (!has_1d) {
          LOG_ERROR("union mode requires --1d-coords-file");
          return 1;
        }
        for (std::size_t nid = 0; nid < node_1d_coords.size(); ++nid) {
          double node_start = node_1d_coords[nid];
          double node_end = node_start + aln_graph.node(nid).sequence.size();
          if (node_end > coord_1d_start && node_start < coord_1d_end) {
            roi_nodes.insert(nid);
          }
        }
      }

      // Add reverse counterparts
      for (std::size_t nid : roi_nodes) {
        roi_with_rev.insert(nid);
        std::size_t orig = piru::index::originalIndex(nid);
        if (piru::index::isReverseNode(nid)) {
          roi_with_rev.insert(piru::index::forwardNodeId(orig));
        } else {
          roi_with_rev.insert(piru::index::reverseNodeId(orig));
        }
      }
    }

    LOG_INFO("  1d=[" + std::to_string(coord_1d_start) + ", " +
             std::to_string(coord_1d_end) + ")" +
             (with_nodes ? " " + mode + " nodes: " + std::to_string(roi_with_rev.size()) : ""));

    /* Write output line */
    // Default: region_label\t1d_start\t1d_end
    // With --with-nodes: region_label\t1d_start\t1d_end\tnode1,node2,...
    ofs << rec.path_name << ":" << rec.start << "-" << rec.end;
    if (has_1d) {
      ofs << "\t" << coord_1d_start << "\t" << coord_1d_end;
    } else {
      ofs << "\t*\t*";
    }
    if (with_nodes) {
      ofs << "\t";
      bool first = true;
      for (std::size_t nid : roi_with_rev) {
        if (!first) ofs << ",";
        ofs << nid;
        first = false;
      }
    }
    ofs << "\n";
  }

  LOG_INFO("Annotation written to " + output_path +
           (with_nodes ? " (" + mode + " nodes)" : " (1D intervals only)"));
  return 0;
}
