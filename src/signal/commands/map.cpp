/**
 * map.cpp
 *
 * CLI handler for `piru map`. Loads a pre-built index, processes reads
 * through the mapping pipeline (seeds -> anchors -> chains -> results),
 * and writes output in PAF/GAF/GAM/JSON format.
 *
 * Related:
 *  - batch_mapper.cpp    (mapping pipeline)
 *  - serialization.cpp   (.pirx loading)
 *
 * SPDX-License-Identifier: MIT
 */

#include "signal/commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "core/index/sort_1d.hpp"
#include "core/io/index/serialization.hpp"
#include "signal/io/reads/read_provider_factory.hpp"
#include "core/io/results/result_writer_factory.hpp"
#include "signal/mapping/batch_mapper.hpp"
#include "core/mapping/path_chainer.hpp"
#include "core/util/logging.hpp"
#include "core/util/timing.hpp"
#include "version.hpp"

namespace {

panomap::signal::SeedExtractorConfig seed_config_from_store(const panomap::index::SeedStore& store) {
  panomap::signal::SeedExtractorConfig cfg;
  cfg.backend = store.extractor_name();
  const auto& params = store.params();
  auto get_u64 = [&params](const std::string& key, std::size_t default_val) -> std::size_t {
    auto it = params.find(key);
    if (it == params.end()) return default_val;
    return std::stoull(it->second);
  };
  cfg.k = get_u64("k", cfg.k);
  cfg.stride = get_u64("stride", cfg.stride);
  cfg.window = get_u64("window", cfg.window);
  auto it_qbits = params.find("qbits");
  if (it_qbits != params.end())
    cfg.qbits = static_cast<std::uint32_t>(std::stoul(it_qbits->second));
  auto it_params = params.find("params");
  if (it_params != params.end()) cfg.params = it_params->second;
  return cfg;
}

}  // namespace

int handle_map(const std::vector<std::string>& args) {
  /* CLI parsing */

  panomap::cli::Parsed parsed;
  panomap::cli::ParseConfig config;
  config.usage = "Usage: piru map [options] --index <file> <reads-path>";
  config.positional_help = {
      "<reads-path>       Input slow5/blow5 file or directory containing reads"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'i', "index", true, "Path to pre-built index file (.pirx)"},
      {'t', "threads", true, "Worker threads (default: 1)"},
      {'p', "profile", false, "Emit timing profile (tree)"},
      {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
      {'\0', "", false, "\nSignal Processing Options:"},
      {'\0', "chunk-size", true, "Signal chunk size in samples (default: 4000, 0 = no chunking)"},
      {'\0', "max-chunks", true, "Max chunks to process per read (default: 10, 0 = unlimited)"},
      {'\0', "sensitivity", true, "Event detection sensitivity (default: 1.0, higher = more sensitive)"},
      {'\0', "diff", true,
       "Diff filter: skip events within diff of last emitted (default: 0.35)"},
      {'\0', "landmark-prominence", true,
       "Landmark: min peak prominence in normalized units (default: 0.5; "
       "must match index-time setting)"},
      // Hidden: individual event params (still work, override sensitivity)
      {'\0', "event-w1", true, ""},
      {'\0', "event-w2", true, ""},
      {'\0', "event-t1", true, ""},
      {'\0', "event-t2", true, ""},
      {'\0', "event-peak", true, ""},
      {'\0', "", false, "\nSeed Lookup Options:"},
      {'\0', "seed-freq-cap", true,
       "Skip high-freq seeds: <N> for absolute, p<F> for percentile (default: p0.99)"},
      {'\0', "max-hits", true,
       "Max seed hits per read before stopping lookup (default: 100000, 0 = unlimited)"},
      {'\0', "", false, "\nChaining Options:"},
      {'\0', "chainer", true,
       "Chainer backend: path-chain (default), sort-chain, pan-chain"},
      {'\0', "1d-coords-file", true,
       "1D coords TSV for sort-chain/pan-chain (from odgi sort --path-sgd-layout)"},
  };
  // Append backend-specific chaining options
  auto chain_opts = panomap::mapping::PathChainerConfig::cli_options();
  config.options.insert(config.options.end(), chain_opts.begin(), chain_opts.end());
  config.options.push_back(
      {'\0', "chain-dd-tolerance", true,
       "SortChainer: dead zone fraction for diagonal deviation (default: 0.0)"});
  config.options.push_back({'\0', "chain-pen-ratio", true,
                            "SortChainer: ratio consistency penalty factor (default: 0.5)"});
  config.options.push_back({'\0', "chain-bidirectional", false,
                            "SortChainer: run DP in both 1D directions, merge chains"});
  config.options.push_back({'\0', "chain-band-1d", true,
                            "PanChainer: 1D band width for candidate selection (default: 5000)"});
  config.options.push_back({'\0', "chain-pen-switch", true,
                            "PanChainer: penalty for switching haplotype path (default: 50)"});
  config.options.push_back({'\0', "", false, "\nMapping Decision Options:"});
  config.options.push_back({'\0', "map-w-threshold", true, "Multi-chain: weighted standout threshold to call mapped (default: 0.45)"});
  config.options.push_back({'\0', "map-sc-min-anchors", true, "Single-chain: min anchors to accept (default: 5)"});
  config.options.push_back({'\0', "map-sc-ratio-lo", true, "Single-chain: min event/ref ratio (default: 0.7)"});
  config.options.push_back({'\0', "map-sc-ratio-hi", true, "Single-chain: max event/ref ratio (default: 1.4)"});
  config.options.push_back({'\0', "no-early-exit", false, "Disable early exit: process all chunks before deciding"});
  config.options.push_back({'\0', "map-fallback-floor", true, "Fallback: absolute score floor (default: 100). Threshold = max(ema_mean, floor)."});
  config.options.push_back({'\0', "no-adaptive", false, "Disable EMA assist, use floor score only"});
  config.options.push_back({'\0', "", false, "\nOutput Options:"});
  config.options.push_back({'o', "output", true, "Output file (.paf or .gaf, default: GAF to stdout)"});
  config.options.push_back({'\0', "secondary", false, "Output secondary alignments (default: primary only)"});
  config.options.push_back({'\0', "", false, "\nDebug Options:"});
  config.options.push_back({'\0', "dump-seed-store", true, "Dump full seed store hash table to TSV file"});
  config.options.push_back({'\0', "no-anchor-merge", false, "Disable anchor merging (for heatmap debugging)"});

  config.on_error = [](const std::string&) { std::cerr << "map: invalid option\n"; };

  if (!panomap::cli::parse_args(args, config, parsed)) {
    panomap::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    panomap::cli::print_help(config, std::cout);
    return 0;
  }

  /* Argument validation */

  // Validate required --index
  if (!parsed.values.count("index")) {
    LOG_ERROR("map: must specify --index <file>");
    panomap::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.positionals.size() != 1) {
    LOG_ERROR("map: missing required <reads-path>");
    panomap::cli::print_help(config, std::cerr);
    return 1;
  }

  // Extract common options
  const bool profile = parsed.values.count("profile") > 0;
  const bool verbose = parsed.values.count("verbose") > 0;
  if (verbose) {
    panomap::logger.set_level(panomap::LogLevel::DEBUG);
  }
  const int num_threads = [&]() {
    auto it = parsed.values.find("threads");
    if (it == parsed.values.end()) return 1;
    try {
      return std::stoi(it->second);
    } catch (...) {
      LOG_WARN("map: invalid --threads value '" + it->second + "', using auto");
      return -1;
    }
  }();

  // Create executor for parallel operations
  auto executor = panomap::concurrency::make_executor(num_threads);
  LOG_DEBUG("Using " + std::to_string(executor->max_concurrency()) + " threads (" +
            executor->backend_name() + ")");

  // Extract output options
  const std::string output_path = parsed.values.count("output") ? parsed.values.at("output") : "";

  // Result writer created after index load (needs graph for GAF path lookups)

  PANOMAP_PROFILE_START(profile, "map");

  /* Index loading */

  PANOMAP_PROFILE_START(profile, "index");

  const std::string index_path = parsed.values.at("index");
  LOG_INFO("loading index from " + index_path);

  if (!panomap::io::index::is_pirx_index(index_path)) {
    LOG_ERROR("map: unsupported index format '" + index_path + "' (expected .pirx file)");
    return 1;
  }

  auto loaded = panomap::io::index::load_index(index_path);

  if (loaded.metadata.mode != panomap::io::index::IndexMode::kSignal) {
    LOG_ERROR(std::string("map: index was built in mode '") +
              panomap::io::index::mode_name(loaded.metadata.mode) +
              "', but piru-signal only loads 'signal' indexes. Use piru-base instead.");
    return 1;
  }

  LOG_INFO("index loaded: " + std::to_string(loaded.graph->nodeCount()) + " nodes");
  LOG_INFO("index metadata: model=" + loaded.metadata.model_name + ", tokenizer=" +
           loaded.metadata.tokenizer + ", seeds=" + loaded.seeds->extractor_name());

  // Extract path lengths from loaded graph BEFORE moving (AdjListGraphStore has path access)
  std::vector<std::size_t> path_lengths;
  {
    const auto& fg = loaded.graph->flat();
    path_lengths.resize(fg.pathCount());
    for (std::uint32_t i = 0; i < fg.pathCount(); ++i) {
      path_lengths[i] = fg.pathLength(i);
    }
  }

  bool primary_only = !parsed.values.count("secondary");
  const auto& flat_graph = loaded.graph->flat();

  auto graph_store = std::move(loaded.graph);
  auto seed_store = std::move(loaded.seeds);
  auto linearization_coords = std::move(loaded.linearization_coords);
  const std::string tokenizer_name = loaded.metadata.tokenizer;
  const std::string pore_model_name = loaded.metadata.model_name;
  const std::size_t loaded_pore_k = loaded.metadata.pore_k;

  PANOMAP_PROFILE_STOP(profile, "index");

  /* Read file discovery */

  const std::string reads_path = parsed.positionals[0];
  std::vector<std::filesystem::path> files;
  std::error_code ec;
  if (std::filesystem::is_directory(reads_path, ec)) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(reads_path)) {
      if (!entry.is_regular_file()) continue;
      const auto ext = entry.path().extension().string();
      const std::string ext_lower = [&]() {
        std::string s = ext;
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return s;
      }();
      if (ext_lower == ".blow5" || ext_lower == ".slow5") {
        files.push_back(entry.path());
      }
    }
  } else {
    files.push_back(reads_path);
  }

  if (files.empty()) {
    LOG_ERROR("map: no supported read files found under '" + reads_path +
              "' (expects .blow5/.slow5)");
    return 1;
  }

  /* Mapper configuration */

  panomap::mapping::BatchMapperConfig map_config;
  map_config.num_threads = num_threads;
  map_config.seed_store = seed_store.get();
  map_config.graph_store = graph_store.get();
  map_config.linearization_coords = &linearization_coords;  // For DP chaining
  map_config.path_lengths = &path_lengths;                  // For anchor bounds checking

  // Load 1D coords: prefer .pirx embedded, fallback to --1d-coords-file
  std::vector<float> node_1d_coords;
  if (!loaded.node_1d_coords.empty()) {
    node_1d_coords = std::move(loaded.node_1d_coords);
    map_config.node_1d_coords = &node_1d_coords;
  } else if (parsed.values.count("1d-coords-file")) {
    std::size_t num_nodes = graph_store->nodeCount();
    node_1d_coords =
        panomap::index::import_1d_coords_odgi(parsed.values.at("1d-coords-file"), num_nodes);
    map_config.node_1d_coords = &node_1d_coords;
  }

  // Load component IDs from index
  std::vector<std::uint32_t> component_ids;
  if (!loaded.component_ids.empty()) {
    component_ids = std::move(loaded.component_ids);
    map_config.component_ids = &component_ids;
  }

  // Create result writer(s).
  //   - output_path empty -> stdout
  //   - output_path is an existing directory OR ends with '/' -> per-file writers
  //     (one output per input file, named <basename>.gaf inside the dir)
  //   - otherwise -> single writer to the given path
  const std::vector<float>* writer_1d = node_1d_coords.empty() ? nullptr : &node_1d_coords;
  const std::vector<std::uint32_t>* writer_cc = component_ids.empty() ? nullptr : &component_ids;
  panomap::io::ResultWriterPtr result_writer;                    // single-file mode
  std::vector<panomap::io::ResultWriterPtr> per_file_writers;    // per-file mode (owning)
  std::unordered_map<std::string, panomap::io::ResultWriter*> per_file_writer_map;

  auto is_per_file_output = [&]() {
    if (output_path.empty()) return false;
    std::error_code ec;
    if (std::filesystem::is_directory(output_path, ec)) return true;
    if (!output_path.empty() && output_path.back() == '/') return true;
    return false;
  }();

  if (output_path.empty()) {
    result_writer =
        panomap::io::make_result_writer_stdout(flat_graph, primary_only, writer_1d, writer_cc);
  } else if (is_per_file_output) {
    // Ensure output dir exists.
    std::error_code ec;
    std::filesystem::create_directories(output_path, ec);
    // Determine output extension from first input file (default .gaf).
    const std::string out_ext = ".gaf";
    for (const auto& f : files) {
      const std::string base = f.stem().string();  // filename without extension
      std::string out_file = output_path;
      if (out_file.back() != '/') out_file += '/';
      out_file += base + out_ext;
      auto writer = panomap::io::make_result_writer(out_file, flat_graph, primary_only,
                                                 writer_1d, writer_cc);
      per_file_writer_map[f.string()] = writer.get();
      per_file_writers.push_back(std::move(writer));
      LOG_INFO("Per-file output: " + f.string() + " -> " + out_file);
    }
  } else {
    result_writer =
        panomap::io::make_result_writer(output_path, flat_graph, primary_only, writer_1d, writer_cc);
    LOG_INFO("Writing results to: " + output_path);
  }

  // Configure tokenizer to match index-time settings.
  // Must mirror IndexPipelineConfig defaults so index and query produce
  // identical tokens.
  if (!tokenizer_name.empty()) {
    map_config.tokenizer_config.backend = tokenizer_name;
  }
  map_config.tokenizer_config.pore_model = pore_model_name;
  map_config.tokenizer_config.fine_min = -2.0f;
  map_config.tokenizer_config.fine_max = 2.0f;
  map_config.tokenizer_config.fine_range = 0.4f;  // must match IndexPipelineConfig::tokenizer_fine_range
  map_config.tokenizer_config.n_bins =
      0;  // must match IndexPipelineConfig::tokenizer_n_bins (0 = use qbits)
  if (parsed.values.count("landmark-prominence")) {
    map_config.tokenizer_config.landmark_min_prominence =
        std::stof(parsed.values.at("landmark-prominence"));
  }
  if (parsed.values.count("diff")) {
    map_config.diff_filter = std::stof(parsed.values.at("diff"));
  }
  LOG_DEBUG("Using tokenizer: " + map_config.tokenizer_config.backend);

  // Configure seed extractor from seed store parameters
  auto index_seed_cfg = seed_config_from_store(*seed_store);
  if (index_seed_cfg.backend.empty()) {
    LOG_ERROR("map: seed store missing extractor name");
    return 1;
  }

  LOG_DEBUG("Using seed extractor settings: backend=" + index_seed_cfg.backend + ", k=" +
            std::to_string(index_seed_cfg.k) + ", stride=" + std::to_string(index_seed_cfg.stride) +
            ", window=" + std::to_string(index_seed_cfg.window) +
            ", qbits=" + std::to_string(index_seed_cfg.qbits));
  map_config.seed_config = index_seed_cfg;

  // Configure chainer
  if (parsed.values.count("chainer")) {
    map_config.chainer_backend = parsed.values.at("chainer");
  }
  map_config.chainer_parsed = parsed;

  // pore_k from index metadata, used for chainer scoring span
  map_config.pore_k = loaded_pore_k;

  // Configure event pipeline
  map_config.event_pipeline_config.pore_model = pore_model_name;
  if (parsed.values.count("sensitivity")) {
    map_config.event_pipeline_config.sensitivity = std::stof(parsed.values.at("sensitivity"));
  }
  // Event detection parameter overrides (take precedence over sensitivity scaling)
  if (parsed.values.count("event-w1")) {
    map_config.event_pipeline_config.override_window_length1 =
        std::stoi(parsed.values.at("event-w1"));
  }
  if (parsed.values.count("event-w2")) {
    map_config.event_pipeline_config.override_window_length2 =
        std::stoi(parsed.values.at("event-w2"));
  }
  if (parsed.values.count("event-t1")) {
    map_config.event_pipeline_config.override_threshold1 = std::stof(parsed.values.at("event-t1"));
  }
  if (parsed.values.count("event-t2")) {
    map_config.event_pipeline_config.override_threshold2 = std::stof(parsed.values.at("event-t2"));
  }
  if (parsed.values.count("event-peak")) {
    map_config.event_pipeline_config.override_peak_height =
        std::stof(parsed.values.at("event-peak"));
  }

  // Chunk size for signal processing
  if (parsed.values.count("chunk-size")) {
    map_config.event_pipeline_config.chunk_size = std::stoull(parsed.values.at("chunk-size"));
  }
  if (parsed.values.count("max-chunks")) {
    map_config.event_pipeline_config.max_chunks = std::stoull(parsed.values.at("max-chunks"));
  }

  // Per-read total hit cap
  if (parsed.values.count("max-hits")) {
    map_config.max_total_hits = std::stoull(parsed.values.at("max-hits"));
  }

  // Map-time frequency cap: set max seed frequency directly
  {
    // Seed frequency cap: p<float> for percentile, <int> for absolute.
    // Default: p0.99 (99th percentile, matching RH2's mid_occ_frac=0.01)
    std::string freq_cap_str = "p0.99";
    if (parsed.values.count("seed-freq-cap")) {
      freq_cap_str = parsed.values.at("seed-freq-cap");
    }

    if (freq_cap_str.size() > 1 && freq_cap_str[0] == 'p') {
      // Percentile mode
      double percentile = std::stod(freq_cap_str.substr(1));
      seed_store->recompute_threshold_from_percentile(percentile);
      LOG_INFO("Seed freq cap: " + freq_cap_str +
               " -> threshold=" + std::to_string(seed_store->frequency_threshold()));
    } else {
      // Absolute mode
      seed_store->set_frequency_threshold(std::stoull(freq_cap_str));
      LOG_INFO("Seed freq cap: " + std::to_string(seed_store->frequency_threshold()));
    }
  }

  LOG_INFO("Seed frequency threshold: " + std::to_string(seed_store->frequency_threshold()) +
           " (seeds with freq > " + std::to_string(seed_store->frequency_threshold()) +
           " will be skipped)");

  map_config.result_writer = result_writer.get();  // null in per-file mode
  map_config.per_file_writers = per_file_writer_map;

  if (parsed.values.count("no-anchor-merge")) {
    map_config.enable_anchor_merge = false;
    LOG_INFO("Anchor merging disabled");
  }

  /* Mapping decision params */
  if (parsed.values.count("map-w-threshold"))
    map_config.map_w_threshold = std::stof(parsed.values.at("map-w-threshold"));
  if (parsed.values.count("map-sc-min-anchors"))
    map_config.map_sc_min_anchors = std::stoull(parsed.values.at("map-sc-min-anchors"));
  if (parsed.values.count("map-sc-ratio-lo"))
    map_config.map_sc_ratio_lo = std::stof(parsed.values.at("map-sc-ratio-lo"));
  if (parsed.values.count("map-sc-ratio-hi"))
    map_config.map_sc_ratio_hi = std::stof(parsed.values.at("map-sc-ratio-hi"));
  if (parsed.values.count("no-early-exit"))
    map_config.no_early_exit = true;
  if (parsed.values.count("map-fallback-floor"))
    map_config.map_fallback_floor_score = std::stod(parsed.values.at("map-fallback-floor"));
  if (parsed.values.count("no-adaptive"))
    map_config.map_fallback_adaptive = false;
  LOG_INFO("Mapping decision: w-threshold=" + std::to_string(map_config.map_w_threshold) +
           " sc-min-anchors=" + std::to_string(map_config.map_sc_min_anchors) +
           " sc-ratio=[" + std::to_string(map_config.map_sc_ratio_lo) + "," +
           std::to_string(map_config.map_sc_ratio_hi) + "]" +
           " fallback-floor=" + std::to_string(map_config.map_fallback_floor_score) +
           " adaptive=" + (map_config.map_fallback_adaptive ? "on" : "off") +
           " early-exit=" + (map_config.no_early_exit ? "off" : "on"));

  /* Read processing */

  PANOMAP_PROFILE_START(profile, "mapping");

  // One BatchMapper per input file. Each file gets its own EMA state so the
  // adaptive fallback isn't corrupted by cross-file contamination (an artifact
  // of eval setups where on-target and off-target reads are in separate files).
  // In a real mixed stream, interleaving naturally produces balanced EMA;
  // sequential per-file processing requires per-file reset to approximate that.
  std::size_t total_reads = 0;
  std::size_t total_batches = 0;
  std::size_t files_processed = 0;

  for (const auto& f : files) {
    auto provider = panomap::io::make_read_provider(f.string());
    if (!provider) {
      LOG_WARN("map: unsupported read format for '" + f.string() + "', skipping");
      continue;
    }
    ++files_processed;
    panomap::mapping::BatchMapper mapper(*provider, map_config, std::cout);
    const auto stats = mapper.process_all();
    total_reads += stats.reads_processed;
    total_batches += stats.batches;
  }

  PANOMAP_PROFILE_STOP(profile, "mapping");

  if (files_processed == 0) {
    LOG_ERROR("map: no readable input files");
    return 1;
  }

  LOG_INFO("map: done. files=" + std::to_string(files_processed) +
           ", batches=" + std::to_string(total_batches) + ", reads=" + std::to_string(total_reads));

  PANOMAP_PROFILE_STOP(profile, "map");
  if (profile) panomap::timing::report(std::cerr);
  return 0;
}
