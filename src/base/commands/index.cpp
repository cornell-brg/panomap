/**
 * index.cpp (base mode)
 *
 * CLI handler for `piru-base index`. Parses arguments, loads a graph,
 * runs the base-mode minimizer indexer, and serializes the result to a
 * .pirx file with mode=base stamped in the header.
 *
 * Related:
 *  - base/index/base_indexer.cpp  (indexing logic)
 *  - core/io/index/serialization.cpp  (.pirx format, mode byte)
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/commands/index.hpp"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "base/index/base_indexer.hpp"
#include "cli/parse.hpp"
#include "core/concurrency/executor.hpp"
#include "core/index/flat_graph.hpp"
#include "core/index/simple_expand.hpp"
#include "core/index/sort_1d.hpp"
#include "core/io/graphs/graph_loader_factory.hpp"
#include "core/io/index/serialization.hpp"
#include "core/util/logging.hpp"
#include "core/util/timing.hpp"
#include "version.hpp"

int handle_base_index(const std::vector<std::string>& args) {
  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru-base index [options] <graph-file>";
  config.positional_help = {"<graph-file>      Graph file to index"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
      {'o', "output", true, "Output index file (default: <graph-file>.pirx)"},
      {'t', "threads", true, "Worker threads (default: 1)"},
      {'p', "profile", false, "Emit timing profile (tree)"},
      {'k', "kmer", true, "Minimizer k-mer size (default: 15)"},
      {'w', "window", true, "Minimizer window size (default: 10)"},
      {'\0', "no-1d-sort", false, "Skip 1D canonical coordinate computation"},
  };
  config.on_error = [](const std::string&) { std::cerr << "index: invalid option\n"; };

  if (!piru::cli::parse_args(args, config, parsed)) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    piru::cli::print_help(config, std::cout);
    return 0;
  }

  const bool profile = parsed.values.count("profile") > 0;
  const bool verbose = parsed.values.count("verbose") > 0;
  if (verbose) {
    piru::logger.set_level(piru::LogLevel::DEBUG);
  }

  const int num_threads = [&]() {
    auto it = parsed.values.find("threads");
    if (it == parsed.values.end()) return 1;
    try {
      return std::stoi(it->second);
    } catch (...) {
      LOG_WARN("index: invalid --threads value '" + it->second + "', using auto");
      return -1;
    }
  }();
  auto executor = piru::concurrency::make_executor(num_threads);
  LOG_DEBUG("Using " + std::to_string(executor->max_concurrency()) + " threads (" +
            executor->backend_name() + ")");

  if (parsed.positionals.empty()) {
    LOG_ERROR("index: missing required <graph-file>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  const std::string graph_path = parsed.positionals.front();

  const std::size_t k =
      parsed.values.count("kmer") ? std::stoul(parsed.values.at("kmer")) : std::size_t{15};
  const std::size_t w =
      parsed.values.count("window") ? std::stoul(parsed.values.at("window")) : std::size_t{10};
  const bool compute_1d_sort = !parsed.values.count("no-1d-sort");

  std::string output_base = std::filesystem::path(graph_path).stem().string();
  if (parsed.values.count("output")) {
    output_base = parsed.values.at("output");
  }
  std::string output_path = output_base;
  if (!output_path.ends_with(".pirx")) {
    output_path += ".pirx";
  }

  LOG_INFO("input: " + graph_path);
  LOG_INFO("seeds: base_minimizer, k=" + std::to_string(k) + ", w=" + std::to_string(w));
  LOG_INFO("output: " + output_path);

  PIRU_PROFILE_START(profile, "index");

  /* 1. Load graph. */

  auto loader = piru::io::make_graph_loader(graph_path);
  if (!loader) {
    LOG_ERROR("index: unsupported graph format for '" + graph_path + "'");
    return 1;
  }
  piru::io::ImportedGraph imported;
  if (!loader->load(imported)) {
    LOG_ERROR("index: failed to read graph file '" + graph_path + "'");
    return 1;
  }
  LOG_INFO("loaded graph: " + std::to_string(imported.nodes.size()) + " nodes, " +
           std::to_string(imported.edges.size()) + " edges, " +
           std::to_string(imported.paths.size()) + " paths (" + loader->get_format_name() + ")");

  /* 2. Expand to directional FlatGraph. */

  auto stage_start = std::chrono::high_resolution_clock::now();
  auto flat_graph = piru::index::simpleExpandFlat(imported);
  imported = piru::io::ImportedGraph{};
  auto stage_elapsed =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - stage_start)
          .count();
  LOG_INFO("[1/3] Transformed graph to directional graph: " +
           std::to_string(flat_graph.nodeCount()) + " nodes, " +
           std::to_string(flat_graph.edgeCount()) + " edges, " +
           std::to_string(flat_graph.pathCount()) + " paths [" + std::to_string(stage_elapsed) +
           "s]");

  /* 3. Base-mode minimizer indexing. */

  stage_start = std::chrono::high_resolution_clock::now();
  piru::base::BaseBucketIndexConfig bi_config;
  bi_config.k = k;
  bi_config.w = w;
  bi_config.executor = executor.get();
  auto bi_result = piru::base::bucketIndexBase(flat_graph, bi_config);
  stage_elapsed =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - stage_start)
          .count();
  LOG_INFO("[2/3] indexed: " + std::to_string(bi_result.seed_store->size()) + " unique seeds [" +
           std::to_string(stage_elapsed) + "s]");

  for (std::size_t i = 0; i < flat_graph.pathCount(); ++i) {
    flat_graph.setPathLength(static_cast<std::uint32_t>(i), bi_result.path_lengths[i]);
  }

  auto graph_store = std::make_unique<piru::index::FlatGraphStore>(std::move(flat_graph));
  auto component_ids = piru::index::compute_components(graph_store->flat());

  std::vector<float> node_1d_coords;
  if (compute_1d_sort) {
    auto sort_start = std::chrono::high_resolution_clock::now();
    piru::index::Sort1DConfig sort_cfg;
    sort_cfg.num_threads = static_cast<std::size_t>(num_threads > 0 ? num_threads : 1);
    node_1d_coords = piru::index::compute_1d_sort(graph_store->flat(), bi_result.linearization_coords,
                                                  bi_result.path_lengths, sort_cfg, component_ids);
    auto sort_elapsed =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - sort_start)
            .count();
    LOG_INFO("[3/3] 1D sort: " + std::to_string(node_1d_coords.size()) + " node positions [" +
             std::to_string(sort_elapsed) + "s]");
  }

  /* 4. Serialize. */

  PIRU_PROFILE_START(profile, "serialize");

  piru::io::index::IndexMetadata metadata;
  metadata.version = PIRU_VERSION;
  metadata.build_timestamp = static_cast<uint64_t>(
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
  metadata.mode = piru::io::index::IndexMode::kBase;
  metadata.model_name = "";        // signal-only field
  metadata.pore_k = 0;             // signal-only field
  metadata.tokenizer = "";         // signal-only field

  piru::io::index::save_index(output_path, *graph_store, *bi_result.seed_store,
                              bi_result.linearization_coords, metadata, node_1d_coords,
                              component_ids);

  LOG_INFO("Index written to " + output_path);

  PIRU_PROFILE_STOP(profile, "serialize");
  PIRU_PROFILE_STOP(profile, "index");
  if (profile) piru::timing::report(std::cerr);

  return 0;
}
