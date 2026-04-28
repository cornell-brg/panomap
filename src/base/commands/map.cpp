/**
 * map.cpp (base mode)
 *
 * CLI handler for `piru-base map`. Loads a base-mode index (.pirx with
 * mode=base), reads basecalled FASTQ, runs BaseMapper, writes GAF.
 *
 * Single-file input (no directory walking yet). Run twice for on/off-
 * target evals; emits separate output files.
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/commands/map.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "base/io/reads/fastq_provider.hpp"
#include "base/mapping/base_mapper.hpp"
#include "cli/parse.hpp"
#include "core/concurrency/executor.hpp"
#include "core/index/sort_1d.hpp"
#include "core/io/index/serialization.hpp"
#include "core/io/results/result_writer_factory.hpp"
#include "core/mapping/path_chainer.hpp"
#include "core/util/logging.hpp"
#include "core/util/timing.hpp"

int handle_base_map(const std::vector<std::string>& args) {
  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru-base map [options] --index <index.pirx> <reads.fastq[.gz]>";
  config.positional_help = {"<reads.fastq[.gz]>  Basecalled FASTQ input (plain or gzip)"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'i', "index", true, "Path to base-mode index (.pirx)"},
      {'t', "threads", true, "Worker threads (default: 1)"},
      {'p', "profile", false, "Emit timing profile"},
      {'v', "verbose", false, "Enable verbose logging"},
      {'o', "output", true, "Output GAF/PAF (default: GAF to stdout)"},
      {'\0', "secondary", false, "Output secondary alignments (default: primary only)"},
      {'\0', "", false, "\nChunked Evaluation Options (AS-style early decision):"},
      {'\0', "chunk-bp", true,
       "Bases per chunk (default: 0 = whole-read; 450 ~= 1 sec @ R9.4 translocation)"},
      {'\0', "max-chunks", true,
       "Max chunks to process per read (default: 10, 0 = unlimited; only used when chunk-bp > 0)"},
      {'\0', "no-early-exit", false, "Disable early exit; process all chunks before deciding"},
      {'\0', "no-adaptive", false, "Disable EMA assist for fallback; use floor only"},
      {'\0', "", false, "\nSeed Lookup Options (minimap2-style adaptive selection):"},
      {'\0', "mid-occ-frac", true,
       "Soft cap percentile: top fraction of seeds treated as 'high-occ' (default: 2e-4)"},
      {'\0', "min-mid-occ", true, "Lower clamp on auto soft cap (default: 10)"},
      {'\0', "max-mid-occ", true, "Upper clamp on auto soft cap (default: 1000000)"},
      {'\0', "mid-occ", true, "Override auto soft cap with absolute value (default: -1 = auto)"},
      {'\0', "max-occ", true, "Hard cap: drop seeds with hits > N (default: 4095)"},
      {'\0', "occ-dist", true,
       "Adaptive window: keep ~1 high-occ seed per N bp (default: 500, 0 = disable)"},
      {'\0', "seed-freq-cap", true,
       "[legacy] Hard cap as <N> or p<F> percentile (overrides --mid-occ if set)"},
      {'\0', "max-hits", true,
       "Max seed hits per read before stopping lookup (default: 100000, 0 = unlimited)"},
      {'\0', "", false, "\nChaining Options:"},
      {'\0', "chainer", true, "Chainer backend: path-chain (default), sort-chain, pan-chain"},
      {'\0', "no-anchor-merge", false, "Disable anchor merging"},
  };
  // Append backend-specific chaining options (path-chain / sort-chain / pan-chain).
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
  config.options.push_back({'\0', "", false, "\nMapping Decision Options:"});
  config.options.push_back({'\0', "map-w-threshold", true,
                            "Multi-chain: weighted standout threshold (default: 0.45)"});
  config.options.push_back({'\0', "map-sc-min-anchors", true,
                            "Single-chain: min anchors to accept (default: 5)"});
  config.options.push_back({'\0', "map-sc-ratio-lo", true,
                            "Single-chain: min query/ref ratio (default: 0.7)"});
  config.options.push_back({'\0', "map-sc-ratio-hi", true,
                            "Single-chain: max query/ref ratio (default: 1.4)"});
  config.options.push_back({'\0', "map-fallback-floor", true,
                            "Fallback absolute score floor (default: 100)"});
  config.on_error = [](const std::string&) { std::cerr << "map: invalid option\n"; };

  if (!piru::cli::parse_args(args, config, parsed)) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    piru::cli::print_help(config, std::cout);
    return 0;
  }

  if (!parsed.values.count("index")) {
    LOG_ERROR("map: must specify --index <file>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.positionals.size() != 1) {
    LOG_ERROR("map: missing required <reads.fastq[.gz]>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }

  const bool profile = parsed.values.count("profile") > 0;
  if (parsed.values.count("verbose")) piru::logger.set_level(piru::LogLevel::DEBUG);

  const int num_threads = [&]() {
    auto it = parsed.values.find("threads");
    if (it == parsed.values.end()) return 1;
    try {
      return std::stoi(it->second);
    } catch (...) {
      return -1;
    }
  }();

  PIRU_PROFILE_START(profile, "map");

  /* 1. Load index. */

  PIRU_PROFILE_START(profile, "index_load");
  const std::string index_path = parsed.values.at("index");
  if (!piru::io::index::is_pirx_index(index_path)) {
    LOG_ERROR("map: not a .pirx file: " + index_path);
    return 1;
  }
  auto loaded = piru::io::index::load_index(index_path);
  if (loaded.metadata.mode != piru::io::index::IndexMode::kBase) {
    LOG_ERROR(std::string("map: index was built in mode '") +
              piru::io::index::mode_name(loaded.metadata.mode) +
              "', but piru-base only loads 'base' indexes. Use piru-signal instead.");
    return 1;
  }
  LOG_INFO("index loaded: " + std::to_string(loaded.graph->nodeCount()) + " nodes, " +
           std::to_string(loaded.seeds->size()) + " unique seeds (" + loaded.seeds->extractor_name() +
           ")");

  std::vector<std::size_t> path_lengths;
  {
    const auto& fg = loaded.graph->flat();
    path_lengths.resize(fg.pathCount());
    for (std::uint32_t i = 0; i < fg.pathCount(); ++i) path_lengths[i] = fg.pathLength(i);
  }

  const auto& flat_graph = loaded.graph->flat();
  auto graph_store = std::move(loaded.graph);
  auto seed_store = std::move(loaded.seeds);
  auto linearization_coords = std::move(loaded.linearization_coords);
  std::vector<float> node_1d_coords = std::move(loaded.node_1d_coords);
  std::vector<std::uint32_t> component_ids = std::move(loaded.component_ids);
  PIRU_PROFILE_STOP(profile, "index_load");

  /* 2. Build mapper config. */

  piru::base::mapping::BaseMapperConfig cfg;
  cfg.num_threads = num_threads;

  // Pull k/w from index seed-store params so seeder hashes the same way the
  // indexer did. Falls back to defaults if not present (older indexes).
  const auto& seed_params = seed_store->params();
  auto get_sz = [&](const std::string& key, std::size_t def) -> std::size_t {
    auto it = seed_params.find(key);
    if (it == seed_params.end()) return def;
    return std::stoull(it->second);
  };
  cfg.seeder.k = get_sz("k", 15);
  cfg.seeder.w = get_sz("window", 10);
  LOG_INFO("seeder: k=" + std::to_string(cfg.seeder.k) + ", w=" + std::to_string(cfg.seeder.w));

  cfg.seed_store = seed_store.get();
  cfg.graph_store = graph_store.get();
  cfg.linearization_coords = &linearization_coords;
  cfg.path_lengths = &path_lengths;
  cfg.node_1d_coords = node_1d_coords.empty() ? nullptr : &node_1d_coords;
  cfg.component_ids = component_ids.empty() ? nullptr : &component_ids;

  if (parsed.values.count("chainer")) cfg.chainer_backend = parsed.values.at("chainer");
  cfg.chainer_parsed = parsed;
  if (parsed.values.count("no-anchor-merge")) cfg.enable_anchor_merge = false;

  if (parsed.values.count("max-hits")) cfg.max_total_hits = std::stoull(parsed.values.at("max-hits"));

  /* Chunked evaluation. Default chunk_bp=0 keeps whole-read behaviour;
   * if user passes --chunk-bp, default --max-chunks to 10 (signal-mode parity). */
  if (parsed.values.count("chunk-bp")) {
    cfg.chunk_bp = std::stoull(parsed.values.at("chunk-bp"));
    cfg.max_chunks = 10;
  }
  if (parsed.values.count("max-chunks")) {
    cfg.max_chunks = std::stoull(parsed.values.at("max-chunks"));
  }
  if (parsed.values.count("no-early-exit")) cfg.no_early_exit = true;
  if (parsed.values.count("no-adaptive")) cfg.map_fallback_adaptive = false;

  /* Soft cap (mid_occ). Three ways to set: --mid-occ <N> absolute,
   * --seed-freq-cap (legacy) absolute or p<F>, or --mid-occ-frac <F>
   * (default: 2e-4 minimap2-style auto). */
  {
    if (parsed.values.count("max-occ"))
      cfg.max_max_occ = std::stoull(parsed.values.at("max-occ"));
    if (parsed.values.count("occ-dist"))
      cfg.occ_dist = std::stoull(parsed.values.at("occ-dist"));
    if (parsed.values.count("min-mid-occ"))
      cfg.min_mid_occ = std::stoull(parsed.values.at("min-mid-occ"));
    if (parsed.values.count("max-mid-occ"))
      cfg.max_mid_occ = std::stoull(parsed.values.at("max-mid-occ"));
    if (parsed.values.count("mid-occ-frac"))
      cfg.mid_occ_frac = std::stof(parsed.values.at("mid-occ-frac"));

    if (parsed.values.count("mid-occ")) {
      cfg.mid_occ_override = std::stoll(parsed.values.at("mid-occ"));
    } else if (parsed.values.count("seed-freq-cap")) {
      // Legacy path: parse N or pF, set threshold directly.
      std::string s = parsed.values.at("seed-freq-cap");
      if (s.size() > 1 && s[0] == 'p') {
        seed_store->recompute_threshold_from_percentile(std::stod(s.substr(1)));
      } else {
        seed_store->set_frequency_threshold(std::stoull(s));
      }
      cfg.mid_occ_override = static_cast<std::int64_t>(seed_store->frequency_threshold());
    } else {
      // Auto: mm2-style mid_occ_frac with min/max bounds.
      seed_store->recompute_threshold_from_top_frac(cfg.mid_occ_frac, cfg.min_mid_occ,
                                                    cfg.max_mid_occ);
    }

    LOG_INFO("seed filter: mid_occ=" + std::to_string(seed_store->frequency_threshold()) +
             " max_occ=" + std::to_string(cfg.max_max_occ) +
             " occ_dist=" + std::to_string(cfg.occ_dist) +
             " (frac=" + std::to_string(cfg.mid_occ_frac) + ")");
  }

  /* Mapping decision overrides. */
  if (parsed.values.count("map-w-threshold"))
    cfg.map_w_threshold = std::stof(parsed.values.at("map-w-threshold"));
  if (parsed.values.count("map-sc-min-anchors"))
    cfg.map_sc_min_anchors = std::stoull(parsed.values.at("map-sc-min-anchors"));
  if (parsed.values.count("map-sc-ratio-lo"))
    cfg.map_sc_ratio_lo = std::stof(parsed.values.at("map-sc-ratio-lo"));
  if (parsed.values.count("map-sc-ratio-hi"))
    cfg.map_sc_ratio_hi = std::stof(parsed.values.at("map-sc-ratio-hi"));
  if (parsed.values.count("map-fallback-floor"))
    cfg.map_fallback_floor_score = std::stod(parsed.values.at("map-fallback-floor"));

  LOG_INFO("decision: w-threshold=" + std::to_string(cfg.map_w_threshold) +
           " sc-min-anchors=" + std::to_string(cfg.map_sc_min_anchors) + " sc-ratio=[" +
           std::to_string(cfg.map_sc_ratio_lo) + "," + std::to_string(cfg.map_sc_ratio_hi) + "]");

  /* 3. Output writer. */

  bool primary_only = !parsed.values.count("secondary");
  const std::string output_path =
      parsed.values.count("output") ? parsed.values.at("output") : std::string{};
  const std::vector<float>* writer_1d = node_1d_coords.empty() ? nullptr : &node_1d_coords;
  const std::vector<std::uint32_t>* writer_cc =
      component_ids.empty() ? nullptr : &component_ids;
  piru::io::ResultWriterPtr result_writer;
  if (output_path.empty()) {
    result_writer =
        piru::io::make_result_writer_stdout(flat_graph, primary_only, writer_1d, writer_cc);
  } else {
    result_writer =
        piru::io::make_result_writer(output_path, flat_graph, primary_only, writer_1d, writer_cc);
    LOG_INFO("Writing results to: " + output_path);
  }
  cfg.result_writer = result_writer.get();

  /* 4. Run mapper on the single FASTQ input. */

  const std::string reads_path = parsed.positionals[0];
  piru::base::io::FastqProvider provider(reads_path);
  if (!provider.is_open()) {
    LOG_ERROR("map: failed to open " + reads_path);
    return 1;
  }
  LOG_INFO("input: " + reads_path);

  PIRU_PROFILE_START(profile, "mapping");
  piru::base::mapping::BaseMapper mapper(provider, cfg, std::cout);
  auto stats = mapper.process_all();
  PIRU_PROFILE_STOP(profile, "mapping");

  LOG_INFO("map: done. reads=" + std::to_string(stats.reads_processed) +
           ", mapped=" + std::to_string(stats.reads_mapped) +
           ", unmapped=" + std::to_string(stats.reads_unmapped));

  PIRU_PROFILE_STOP(profile, "map");
  if (profile) piru::timing::report(std::cerr);
  return 0;
}
