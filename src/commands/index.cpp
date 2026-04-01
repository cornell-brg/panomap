/**
 * index.cpp
 *
 * CLI handler for `piru index`. Parses arguments, loads a graph and
 * pore model, runs the indexing pipeline, and serializes the result
 * to a .pirx file.
 *
 * Related:
 *  - index_pipeline.cpp  (indexing logic)
 *  - serialization.cpp   (.pirx format)
 *
 * SPDX-License-Identifier: MIT
 */

#include "commands/index.hpp"

#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/index_pipeline.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/index/serialization.hpp"
#include "io/models/model_factory.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

namespace {

piru::io::ModelPtr load_model_or_file(const std::string& model_arg) {
  if (auto built = piru::io::load_builtin_model(model_arg)) {
    return built;
  }
  if (std::filesystem::exists(model_arg)) {
    return piru::io::load_model_from_file(model_arg);
  }
  LOG_ERROR("Unknown pore model: " + model_arg + " (not a built-in name or file)");
  return nullptr;
}

}  // namespace

int handle_index(const std::vector<std::string>& args) {
  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru index [options] <graph-file>";
  config.positional_help = {"<graph-file>      Graph file to index"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
      {'m', "model", true,
       "Pore model (builtin name: r9.4/r10.4 or model file path, default: r10.4)"},
      {'o', "output", true, "Output index file (default: <graph-file>.pirx)"},
      {'t', "threads", true, "Worker threads"},
      {'p', "profile", false, "Emit timing profile (tree)"},
      {'\0', "", false, "\nSeed Generation Options:"},
      {'\0', "seed-type", true, "Seed extractor type: kmer, minimizer (default)"},
      {'\0', "seed-k", true, "Seed k-mer size (default: 6)"},
      {'\0', "minimizer-window", true,
       "Minimizer window size (default: 5, only with --seed-type minimizer)"},
      {'\0', "seed-stride", true, "Seed stride (default: 1)"},
      {'\0', "seed-freq-cutoff", true, "Seed frequency filter percentile (0.0-1.0, default: 0.9)"},
      {'\0', "seed-freq-cap", true,
       "Subsample cap percentile for filtered seeds (0.0-1.0, default: 0.25)"},
      {'\0', "diff", true,
       "Event diff filter: skip events within diff of last emitted (default: 0, RH2: 0.35)"},
      {'\0', "", false, "\nIndexer Options:"},
      {'\0', "no-1d-sort", false, "Skip 1D canonical coordinate computation"},
      {'\0', "1d-coords-file", true,
       "Import pre-computed 1D coords from TSV (overrides built-in PG-SGD)"},
      {'\0', "", false, "\nDebug Options:"},
      {'\0', "dump-1d-coords", true,
       "Dump 1D sort coordinates to TSV file"},
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

  /* Thread setup */

  const int num_threads = [&]() {
    auto it = parsed.values.find("threads");
    if (it == parsed.values.end()) return -1;
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

  /* Input validation */

  const std::string model_arg = parsed.values.count("model") ? parsed.values.at("model") : "r10.4";

  if (parsed.positionals.empty()) {
    LOG_ERROR("index: missing required <graph-file>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  const std::string graph_path = parsed.positionals.front();

  auto model = load_model_or_file(model_arg);
  if (!model) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }

  PIRU_PROFILE_START(profile, "index");

  auto loader = piru::io::make_graph_loader(graph_path);
  if (!loader) {
    LOG_ERROR("index: unsupported graph format for '" + graph_path + "'");
    piru::cli::print_help(config, std::cerr);
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

  /* Parameter parsing */

  const std::size_t pore_k = model->k();

  // Use defaults from IndexPipelineConfig (single source of truth)
  piru::index::IndexPipelineConfig defaults;
  const std::string seed_type =
      parsed.values.count("seed-type") ? parsed.values.at("seed-type") : defaults.seed_type;
  const std::size_t seed_k =
      parsed.values.count("seed-k") ? std::stoul(parsed.values.at("seed-k")) : defaults.seed_k;
  const std::size_t minimizer_window = parsed.values.count("minimizer-window")
                                           ? std::stoul(parsed.values.at("minimizer-window"))
                                           : defaults.minimizer_window;
  const std::size_t seed_stride = parsed.values.count("seed-stride")
                                      ? std::stoul(parsed.values.at("seed-stride"))
                                      : defaults.seed_stride;
  const double seed_freq_cutoff = parsed.values.count("seed-freq-cutoff")
                                      ? std::stod(parsed.values.at("seed-freq-cutoff"))
                                      : defaults.seed_freq_cutoff;
  const double seed_freq_cap = parsed.values.count("seed-freq-cap")
                                   ? std::stod(parsed.values.at("seed-freq-cap"))
                                   : defaults.seed_freq_cap;
  const float fuzzy_diff = parsed.values.count("diff")
                               ? std::stof(parsed.values.at("diff"))
                               : defaults.fuzzy_diff;

  std::string output_base = std::filesystem::path(graph_path).stem().string();
  if (parsed.values.count("output")) {
    output_base = parsed.values.at("output");
  }

  LOG_INFO("input: " + graph_path);
  LOG_INFO("model: " + model->name() + " (k=" + std::to_string(pore_k) + ")");
  LOG_INFO("seeds: type=" + seed_type + ", k=" + std::to_string(seed_k) + ", window=" +
           std::to_string(minimizer_window) + ", stride=" + std::to_string(seed_stride) +
           ", freq_cutoff=" + std::to_string(seed_freq_cutoff) +
           ", freq_cap=" + std::to_string(seed_freq_cap));
  LOG_INFO("output: " + output_base);

  /* Run indexing pipeline */

  piru::index::IndexPipelineConfig index_config;
  index_config.seed_type = seed_type;
  index_config.seed_k = seed_k;
  index_config.minimizer_window = minimizer_window;
  index_config.seed_stride = seed_stride;
  index_config.seed_freq_cutoff = seed_freq_cutoff;
  index_config.seed_freq_cap = seed_freq_cap;
  index_config.fuzzy_diff = fuzzy_diff;
  index_config.fuzzy_quantizer = "rh2";
  if (parsed.values.count("no-1d-sort")) {
    index_config.compute_1d_sort = false;
  }

  index_config.executor = executor.get();

  auto result = piru::index::run_index_pipeline(imported, *model, index_config);

  // Import pre-computed 1D coords (overrides --compute-1d-sort)
  if (parsed.values.count("1d-coords-file")) {
    std::size_t num_nodes = result.graph_store->nodeCount();
    result.node_1d_coords = piru::index::import_1d_coords_odgi(
        parsed.values.at("1d-coords-file"), num_nodes);
  }

  // Dump 1D coordinates if requested
  if (parsed.values.count("dump-1d-coords") && !result.node_1d_coords.empty()) {
    auto* adj_store = dynamic_cast<piru::index::AdjListGraphStore*>(result.graph_store.get());
    if (adj_store) {
      piru::index::dump_1d_coords_tsv(parsed.values.at("dump-1d-coords"),
                                       result.node_1d_coords,
                                       adj_store->flat());
    }
  }

  /* Serialize index */

  PIRU_PROFILE_START(profile, "serialize");

  std::string output_path = output_base;
  if (!output_path.ends_with(".pirx")) {
    output_path += ".pirx";
  }

  piru::io::index::IndexMetadata metadata;
  metadata.model_name = result.model_name;
  metadata.pore_k = result.pore_k;
  metadata.fuzzy_quantizer = result.fuzzy_quantizer;
  piru::io::index::save_index(output_path, *result.graph_store, *result.seed_store,
                              result.linearization_coords, metadata,
                              result.node_1d_coords);

  LOG_INFO("Index written to " + output_path);

  PIRU_PROFILE_STOP(profile, "serialize");
  PIRU_PROFILE_STOP(profile, "index");
  if (profile) piru::timing::report(std::cerr);

  return 0;
}
