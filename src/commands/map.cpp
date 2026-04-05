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

#include "commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/sort_1d.hpp"
#include "io/index/serialization.hpp"
#include "io/reads/read_provider_factory.hpp"
#include "io/regions/pira_parser.hpp"
#include "io/results/result_writer_factory.hpp"
#include "mapping/batch_mapper.hpp"
#include "mapping/path_chainer.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

namespace {

piru::signal::SeedExtractorConfig seed_config_from_store(const piru::index::SeedStore& store) {
  piru::signal::SeedExtractorConfig cfg;
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

  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru map [options] --index <file> <reads-path>";
  config.positional_help = {
      "<reads-path>       Input slow5/blow5 file or directory containing reads"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'i', "index", true, "Path to pre-built index file (.pirx)"},
      {'t', "threads", true, "Worker threads (-1 = auto)"},
      {'p', "profile", false, "Emit timing profile (tree)"},
      {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
      {'\0', "", false, "\nMapping Options:"},
      {'\0', "seed-freq-cap", true,
       "Skip high-freq seeds: <N> for absolute, p<F> for percentile (default: p0.99)"},
      {'\0', "max-hits", true,
       "Max seed hits per read before stopping lookup (default: 100000, 0 = unlimited)"},
      {'\0', "diff", true,
       "Event diff filter: skip events within diff of last emitted (default: 0, RH2: 0.35)"},
      {'\0', "chunk-size", true, "Signal chunk size in samples (default: 4000, 0 = no chunking)"},
      {'\0', "max-chunks", true, "Max chunks to process per read (default: 10, 0 = unlimited)"},
      {'\0', "chainer", true,
       "Chainer backend: path-chain (default), graph-chain, sort-chain, pan-chain"},
      {'\0', "1d-coords-file", true,
       "1D coords TSV for sort-chain/pan-chain (from odgi sort --path-sgd-layout)"},
      {'\0', "", false, "\nSignal Processing Options:"},
      {'\0', "event-pipeline", true,
       "Event pipeline backend: rawhash (default), scrappie, passthrough"},
      {'\0', "event-w1", true, "Event detection short window length (default: 3)"},
      {'\0', "event-w2", true, "Event detection long window length (default: backend-specific)"},
      {'\0', "event-t1", true, "Event detection threshold1 (default: backend-specific)"},
      {'\0', "event-t2", true, "Event detection threshold2 (default: backend-specific)"},
      {'\0', "event-peak", true, "Event detection peak height (default: backend-specific)"},
      {'\0', "", false, "\nDebug Options:"},
      {'\0', "dump-seed-store", true, "Dump full seed store hash table to TSV file"},
      {'\0', "no-anchor-merge", false, "Disable anchor merging (for heatmap debugging)"},
      {'\0', "", false, "\nClassification Options:"},
      {'\0', "roi", true, "ROI annotation file (.pira from piru annotate)"},
      {'\0', "mode", true, "Classification mode: enrich or deplete (requires --roi)"},
      {'\0', "chain-target", true,
       "ROI-only chaining: filter anchors to ROI nodes, classify by score >= threshold"},
      {'\0', "chain-genome", true,
       "Whole-genome chaining: classify by ROI overlap fraction >= threshold (default: 0.5)"},
      {'\0', "", false, "\nOutput Options:"},
      {'o', "output", true, "Output file path (format auto-detected from extension: .paf, .gaf)"},
      {'\0', "output-format", true, "Override output format (paf, gaf)"},
      {'\0', "secondary", false, "Output secondary alignments (default: primary only)"},
      {'\0', "min-secondary-ratio", true,
       "Min chain score ratio vs primary for secondaries (default: 0.4)"},
      {'\0', "map-min-anchors", true, "Noise floor: min anchors in primary chain (default: 3)"},
      {'\0', "map-min-score", true, "Noise floor: min primary chain score (default: 30)"},
      {'\0', "map-standout-ratio", true, "Fraction of mapq from standout vs score (default: 0.17)"},
      {'\0', "map-min-mapq-exit", true, "Min mapq to call mapped and early exit (default: 12)"},
  };
  // Append backend-specific CLI options
  auto chain_opts = piru::mapping::PathChainerConfig::cli_options();
  config.options.insert(config.options.end(), chain_opts.begin(), chain_opts.end());
  config.options.push_back(
      {'\0', "chain-dd-tolerance", true,
       "SortChainer: dead zone fraction for diagonal deviation (default: 0.0)"});
  config.options.push_back({'\0', "chain-pen-ratio", true,
                            "SortChainer: ratio consistency penalty factor (default: 0.5)"});
  config.options.push_back({'\0', "chain-band-1d", true,
                            "PanChainer: 1D band width for candidate selection (default: 5000)"});
  config.options.push_back({'\0', "chain-pen-switch", true,
                            "PanChainer: penalty for switching haplotype path (default: 50)"});

  config.on_error = [](const std::string&) { std::cerr << "map: invalid option\n"; };

  if (!piru::cli::parse_args(args, config, parsed)) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    piru::cli::print_help(config, std::cout);
    return 0;
  }

  /* Argument validation */

  // Validate required --index
  if (!parsed.values.count("index")) {
    LOG_ERROR("map: must specify --index <file>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.positionals.size() != 1) {
    LOG_ERROR("map: missing required <reads-path>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }

  // Extract common options
  const bool profile = parsed.values.count("profile") > 0;
  const bool verbose = parsed.values.count("verbose") > 0;
  if (verbose) {
    piru::logger.set_level(piru::LogLevel::DEBUG);
  }
  const int num_threads = [&]() {
    auto it = parsed.values.find("threads");
    if (it == parsed.values.end()) return -1;
    try {
      return std::stoi(it->second);
    } catch (...) {
      LOG_WARN("map: invalid --threads value '" + it->second + "', using auto");
      return -1;
    }
  }();

  // Classification options
  const std::string roi_path = parsed.values.count("roi") ? parsed.values.at("roi") : "";
  const std::string classify_mode = parsed.values.count("mode") ? parsed.values.at("mode") : "";
  const bool has_chain_target = parsed.values.count("chain-target") > 0;
  const bool has_chain_genome = parsed.values.count("chain-genome") > 0;

  if (!roi_path.empty() && classify_mode.empty()) {
    LOG_ERROR("map: --roi requires --mode (enrich or deplete)");
    return 1;
  }
  if (!classify_mode.empty() && roi_path.empty()) {
    LOG_ERROR("map: --mode requires --roi <file>");
    return 1;
  }
  if (!classify_mode.empty() && classify_mode != "enrich" && classify_mode != "deplete") {
    LOG_ERROR("map: --mode must be 'enrich' or 'deplete', got '" + classify_mode + "'");
    return 1;
  }
  if (has_chain_target && has_chain_genome) {
    LOG_ERROR("map: --chain-target and --chain-genome are mutually exclusive");
    return 1;
  }
  if ((has_chain_target || has_chain_genome) && roi_path.empty()) {
    LOG_ERROR("map: --chain-target/--chain-genome require --roi <file>");
    return 1;
  }
  if (!roi_path.empty() && !has_chain_target && !has_chain_genome) {
    LOG_ERROR("map: --roi requires --chain-target <score> or --chain-genome <overlap>");
    return 1;
  }

  // Parse classification thresholds
  const bool roi_filter_anchors = has_chain_target;
  const double roi_score_threshold = [&]() {
    if (!has_chain_target) return 0.0;
    try {
      return std::stod(parsed.values.at("chain-target"));
    } catch (...) {
      LOG_WARN("map: invalid --chain-target value, using 30");
      return 30.0;
    }
  }();
  const double roi_overlap_threshold = [&]() {
    if (!has_chain_genome) return 0.5;  // default for implicit genome mode
    try {
      return std::stod(parsed.values.at("chain-genome"));
    } catch (...) {
      LOG_WARN("map: invalid --chain-genome value, using 0.5");
      return 0.5;
    }
  }();

  // Load ROI annotation if provided
  std::unordered_set<std::size_t> roi_nodes;
  if (!roi_path.empty()) {
    roi_nodes = piru::io::parse_pira(roi_path);
  }

  // Create executor for parallel operations (indexing and mapping)
  auto executor = piru::concurrency::make_executor(num_threads);
  LOG_DEBUG("Using " + std::to_string(executor->max_concurrency()) + " threads (" +
            executor->backend_name() + ")");

  // Extract output options
  const std::string output_path = parsed.values.count("output") ? parsed.values.at("output") : "";
  const std::string output_format =
      parsed.values.count("output-format") ? parsed.values.at("output-format") : "";

  // Result writer created after index load (needs graph for GAF path lookups)

  PIRU_PROFILE_START(profile, "map");

  /* Index loading */

  PIRU_PROFILE_START(profile, "index");

  const std::string index_path = parsed.values.at("index");
  LOG_INFO("loading index from " + index_path);

  if (!piru::io::index::is_pirx_index(index_path)) {
    LOG_ERROR("map: unsupported index format '" + index_path + "' (expected .pirx file)");
    return 1;
  }

  auto loaded = piru::io::index::load_index(index_path);

  LOG_INFO("index loaded: " + std::to_string(loaded.graph->nodeCount()) + " nodes");
  LOG_INFO("index metadata: model=" + loaded.metadata.model_name + ", fuzzy=" +
           loaded.metadata.fuzzy_quantizer + ", seeds=" + loaded.seeds->extractor_name());

  // Extract path lengths from loaded graph BEFORE moving (AdjListGraphStore has path access)
  std::vector<std::size_t> path_lengths;
  {
    const auto& fg = loaded.graph->flat();
    path_lengths.resize(fg.pathCount());
    for (std::uint32_t i = 0; i < fg.pathCount(); ++i) {
      path_lengths[i] = fg.pathLength(i);
    }
  }

  // Create result writer (needs graph for GAF path lookups)
  bool primary_only = !parsed.values.count("secondary");

  piru::io::ResultWriterPtr result_writer;
  if (!output_path.empty()) {
    const auto& flat_graph = loaded.graph->flat();
    if (!output_format.empty()) {
      result_writer = piru::io::make_result_writer(output_path, output_format, flat_graph, primary_only);
    } else {
      result_writer = piru::io::make_result_writer(output_path, flat_graph, primary_only);
    }
    if (!result_writer) {
      LOG_ERROR("map: failed to create output writer for '" + output_path + "'");
      return 1;
    }
    LOG_INFO("Writing results to: " + output_path);
  }

  auto graph_store = std::move(loaded.graph);
  auto seed_store = std::move(loaded.seeds);
  auto linearization_coords = std::move(loaded.linearization_coords);
  const std::string fuzzy_quantizer_name = loaded.metadata.fuzzy_quantizer;
  const std::string pore_model_name = loaded.metadata.model_name;
  const std::size_t loaded_pore_k = loaded.metadata.pore_k;

  PIRU_PROFILE_STOP(profile, "index");

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

  piru::mapping::BatchMapperConfig map_config;
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
        piru::index::import_1d_coords_odgi(parsed.values.at("1d-coords-file"), num_nodes);
    map_config.node_1d_coords = &node_1d_coords;
  }

  // Configure fuzzy quantizer to match index-time settings.
  // Must mirror IndexPipelineConfig defaults so index and query produce
  // identical tokens. The factory only overrides for "piru" backend;
  // "rh2" backend uses these values directly.
  if (!fuzzy_quantizer_name.empty()) {
    map_config.fuzzy_config.backend = fuzzy_quantizer_name;
  }
  map_config.fuzzy_config.pore_model = pore_model_name;
  map_config.fuzzy_config.fine_min = -2.0f;
  map_config.fuzzy_config.fine_max = 2.0f;
  map_config.fuzzy_config.fine_range = 0.4f;  // must match IndexPipelineConfig::fuzzy_fine_range
  map_config.fuzzy_config.n_bins =
      0;  // must match IndexPipelineConfig::fuzzy_n_bins (0 = use qbits)
  if (parsed.values.count("diff")) {
    map_config.fuzzy_config.diff = std::stof(parsed.values.at("diff"));
  }
  LOG_DEBUG("Using fuzzy quantizer: " + map_config.fuzzy_config.backend);

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

  // ROI classification config
  if (!roi_nodes.empty()) {
    map_config.roi_nodes = &roi_nodes;
    map_config.classify_mode = classify_mode;
    map_config.roi_filter_anchors = roi_filter_anchors;
    map_config.roi_score_threshold = roi_score_threshold;
    map_config.roi_overlap_threshold = roi_overlap_threshold;
    // Force top-1 chain only -- we only need the best for classification
    map_config.chainer_parsed.values["chain-max-chains"] = "1";
    if (roi_filter_anchors) {
      LOG_INFO("ROI classification: mode=" + classify_mode +
               ", chain-target score>=" + std::to_string(roi_score_threshold) + ", " +
               std::to_string(roi_nodes.size()) + " ROI nodes");
    } else {
      LOG_INFO("ROI classification: mode=" + classify_mode +
               ", chain-genome overlap>=" + std::to_string(roi_overlap_threshold) + ", " +
               std::to_string(roi_nodes.size()) + " ROI nodes");
    }
  }

  // pore_k from index metadata, used for chainer scoring span
  map_config.pore_k = loaded_pore_k;

  // Configure event pipeline (unified event detection + normalization)
  map_config.event_pipeline_config.pore_model = pore_model_name;
  if (parsed.values.count("event-pipeline")) {
    map_config.event_pipeline_config.backend = parsed.values.at("event-pipeline");
  }
  // Event detection parameter overrides (take precedence over backend defaults)
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

  // Configure result writer if specified
  if (result_writer) {
    map_config.result_writer = result_writer.get();
  }

  if (parsed.values.count("no-anchor-merge")) {
    map_config.enable_anchor_merge = false;
    LOG_INFO("Anchor merging disabled");
  }

  // min-secondary-ratio is now handled by the GafWriter directly
  // TODO: pass this to the writer config if needed

  /* Mapping decision params */
  if (parsed.values.count("map-min-anchors"))
    map_config.map_min_anchors = std::stoull(parsed.values.at("map-min-anchors"));
  if (parsed.values.count("map-min-score"))
    map_config.map_min_score = std::stod(parsed.values.at("map-min-score"));
  if (parsed.values.count("map-standout-ratio"))
    map_config.map_standout_ratio = std::stof(parsed.values.at("map-standout-ratio"));
  if (parsed.values.count("map-min-mapq-exit"))
    map_config.map_min_mapq_exit = std::stoi(parsed.values.at("map-min-mapq-exit"));
  LOG_INFO("Mapping decision: min-anchors=" + std::to_string(map_config.map_min_anchors) +
           " min-score=" + std::to_string(map_config.map_min_score) +
           " standout-ratio=" + std::to_string(map_config.map_standout_ratio) +
           " min-mapq-exit=" + std::to_string(map_config.map_min_mapq_exit));

  /* Read processing */

  PIRU_PROFILE_START(profile, "mapping");

  std::size_t total_reads = 0;
  std::size_t total_batches = 0;
  std::size_t files_processed = 0;

  for (const auto& f : files) {
    auto provider = piru::io::make_read_provider(f.string());
    if (!provider) {
      LOG_WARN("map: unsupported read format for '" + f.string() + "', skipping");
      continue;
    }
    ++files_processed;
    piru::mapping::BatchMapper mapper(*provider, map_config, std::cout);
    const auto stats = mapper.process_all();
    total_reads += stats.reads_processed;
    total_batches += stats.batches;
  }

  PIRU_PROFILE_STOP(profile, "mapping");

  if (files_processed == 0) {
    LOG_ERROR("map: no readable input files");
    return 1;
  }

  LOG_INFO("map: done. files=" + std::to_string(files_processed) +
           ", batches=" + std::to_string(total_batches) + ", reads=" + std::to_string(total_reads));

  PIRU_PROFILE_STOP(profile, "map");
  if (profile) piru::timing::report(std::cerr);
  return 0;
}
