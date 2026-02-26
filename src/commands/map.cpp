#include "commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/index_pipeline.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/index/simple_serialization.hpp"
#include "io/models/model_factory.hpp"
#include "io/reads/read_provider_factory.hpp"
#include "io/results/result_writer_factory.hpp"
#include "mapping/batch_mapper.hpp"
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

piru::signal::SeedExtractorConfig seed_config_from_store(const piru::index::HashSeedStore& store) {
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
    if (it_qbits != params.end()) cfg.qbits = static_cast<std::uint32_t>(std::stoul(it_qbits->second));
    auto it_params = params.find("params");
    if (it_params != params.end()) cfg.params = it_params->second;
    return cfg;
}

}  // namespace

int handle_map(const std::vector<std::string>& args) {
    // =========================================================================
    // CLI PARSING
    // =========================================================================

    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru map [options] (--index <dir> | --graph <file> --model <model>) <reads-path>";
    config.positional_help = {"<reads-path>       Input slow5/blow5 file or directory containing reads"};
    config.options = {
        {'h', "help", false, "Show help"},
        {'i', "index", true, "Path to pre-built index directory (XOR with --graph)"},
        {'g', "graph", true, "Graph file for in-memory indexing (XOR with --index)"},
        {'m', "model", true, "Pore model (required with --graph; r9.4/r10.4 or file path)"},
        {'t', "threads", true, "Worker threads (-1 = auto)"},
        {'p', "profile", false, "Emit timing profile (tree)"},
        {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
        {'\0', "", false, "\nMapping Options:"},
        {'\0', "seed-freq-cap", true, "Skip seeds above this frequency percentile at lookup (0.0-1.0, default: none)"},
        {'\0', "clusterer", true, "Clusterer backend: dp-chain (default), fse, probe"},
        {'\0', "chain-max-dist", true, "DP chain: max query/ref distance between anchors (default: 500)"},
        {'\0', "chain-max-diag-dev", true, "DP chain: max diagonal deviation (default: 500)"},
        {'\0', "chain-gap-penalty", true, "DP chain: gap penalty factor (default: 0.02)"},
        {'\0', "chain-diag-penalty", true, "DP chain: diagonal penalty factor (default: 0.05)"},
        {'\0', "chain-overlap-penalty", true, "DP chain: overlap penalty factor (default: 0.4)"},
        {'\0', "chain-anchor-weight", true, "DP chain: anchor weight (default: 1.5)"},
        {'\0', "chain-min-score", true, "DP chain: minimum chain score (default: 12)"},
        {'\0', "chain-max-chains", true, "DP chain: max chains to extract (default: 10)"},
        {'\0', "chain-max-skip", true, "DP chain: stop after N consecutive failed attempts (default: 25)"},
        {'\0', "chain-merge", true, "DP chain: merge overlapping chains (default: true)"},
        {'\0', "", false, "\nSignal Processing Options (only with --graph):"},
        {'\0', "indexer-backend", true, "Indexer backend: node-first, path-walk (default)"},
        {'\0', "event-pipeline", true, "Event pipeline backend: rawhash (default), scrappie, passthrough"},
        {'\0', "event-w1", true, "Event detection short window length (default: 3)"},
        {'\0', "event-w2", true, "Event detection long window length (default: backend-specific)"},
        {'\0', "event-t1", true, "Event detection threshold1 (default: backend-specific)"},
        {'\0', "event-t2", true, "Event detection threshold2 (default: backend-specific)"},
        {'\0', "event-peak", true, "Event detection peak height (default: backend-specific)"},
        {'\0', "fuzzy-backend", true, "Fuzzy quantizer backend: piru (default), rh2"},
        {'\0', "fuzzy-fine-min", true, "Fuzzy quantizer fine region min (default: -2.0)"},
        {'\0', "fuzzy-fine-max", true, "Fuzzy quantizer fine region max (default: 2.0)"},
        {'\0', "fuzzy-fine-range", true, "Fuzzy quantizer fine range (default: 0.9 R9, 0.8 R10)"},
        {'\0', "fuzzy-bins", true, "Fuzzy quantizer bin count (default: 10)"},
        {'\0', "seed-backend", true, "Seed extractor backend: kmer (default), minimizer"},
        {'\0', "seed-k", true, "Seed extractor k-mer size (default: 6)"},
        {'\0', "seed-w", true, "Minimizer window size (default: 5, only with --seed-backend minimizer)"},
        {'\0', "seed-stride", true, "Seed extractor stride (default: 1)"},
        {'\0', "seed-filter", true, "Index-time seed frequency filter percentile (0.0-1.0, default: 1.0)"},
        {'\0', "seed-subsample", true, "Subsample cap percentile for filtered seeds (0.0-1.0, default: 0.55)"},
        {'\0', "seed-mode", true, "Seeding mode: node, path (default)"},
        {'\0', "", false, "\nDebug Options:"},
        {'\0', "dump-anchors", true, "Dump anchors to directory (one file per read)"},
        {'\0', "dump-chains", true, "Dump chains to directory (one file per read)"},
        {'\0', "dump-hit-stats", true, "Dump seed hit statistics to directory (one file per read)"},
        {'\0', "dump-path-chains", true, "Dump best chain per path to directory (diagnostic)"},
        {'\0', "dump-seed-store", true, "Dump full seed store hash table to TSV file"},
        {'\0', "dump-read-seeds", true, "Dump all read seeds (including no-hit) to directory"},
        {'\0', "dump-anchor-detail", true, "Dump per-anchor detail to directory (use with --dump-anchor-reads)"},
        {'\0', "dump-anchor-reads", true, "Comma-separated read ID substrings to dump anchors for"},
        {'\0', "no-anchor-merge", false, "Disable anchor merging (for heatmap debugging)"},
        {'\0', "", false, "\nOutput Options:"},
        {'o', "output", true, "Output file path (format auto-detected from extension: .paf, .gaf, .gam, .json)"},
        {'\0', "output-format", true, "Override output format (paf, gaf, gam, json)"},
        {'\0', "min-secondary-ratio", true, "Min chain score ratio vs primary for secondaries (default: 0.4)"},
    };
    config.on_error = [](const std::string&) { std::cerr << "map: invalid option\n"; };

    if (!piru::cli::parse_args(args, config, parsed)) {
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.values.count("help")) {
        piru::cli::print_help(config, std::cout);
        return 0;
    }

    // =========================================================================
    // ARGUMENT VALIDATION
    // =========================================================================

    // Validate mutually exclusive options: --index XOR --graph
    const bool has_index = parsed.values.count("index") > 0;
    const bool has_graph = parsed.values.count("graph") > 0;

    if (!has_index && !has_graph) {
        LOG_ERROR("map: must specify either --index <dir> or --graph <file>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (has_index && has_graph) {
        LOG_ERROR("map: cannot specify both --index and --graph (choose one)");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (has_graph && !parsed.values.count("model")) {
        LOG_ERROR("map: --graph requires --model <model>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.positionals.size() != 1) {
        LOG_ERROR("map: missing required <reads-path>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }

    // Warn about flags that are ignored with --index (index-time-only params)
    if (has_index) {
        std::vector<std::string> graph_only_flags = {
            "seed-backend", "seed-k", "seed-w", "seed-stride", "seed-filter", "seed-subsample", "seed-mode",
            "indexer-backend",
            "event-pipeline", "event-w1", "event-w2", "event-t1", "event-t2", "event-peak",
            "fuzzy-backend", "fuzzy-fine-min", "fuzzy-fine-max", "fuzzy-fine-range", "fuzzy-bins",
            "model",
        };
        std::vector<std::string> found;
        for (const auto& flag : graph_only_flags) {
            if (parsed.values.count(flag)) found.push_back("--" + flag);
        }
        if (!found.empty()) {
            std::string list;
            for (std::size_t i = 0; i < found.size(); ++i) {
                if (i > 0) list += ", ";
                list += found[i];
            }
            LOG_WARN("Ignoring index-time flags with --index (these only apply with --graph): " + list);
        }
    }

    // Warn about --seed-w without --seed-backend minimizer
    if (parsed.values.count("seed-w") && has_graph) {
        const std::string backend = parsed.values.count("seed-backend")
                                        ? parsed.values.at("seed-backend")
                                        : "kmer";
        if (backend != "minimizer") {
            LOG_WARN("--seed-w is ignored without --seed-backend minimizer");
        }
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

    // Create executor for parallel operations (indexing and mapping)
    auto executor = piru::concurrency::make_executor(num_threads);
    LOG_DEBUG("Using " + std::to_string(executor->max_concurrency()) + " threads (" +
              executor->backend_name() + ")");

    // Extract output options
    const std::string output_path = parsed.values.count("output")
                                        ? parsed.values.at("output")
                                        : "";
    const std::string output_format = parsed.values.count("output-format")
                                          ? parsed.values.at("output-format")
                                          : "";

    // Create result writer if output specified
    piru::io::ResultWriterPtr result_writer;
    if (!output_path.empty()) {
        if (!output_format.empty()) {
            // Use explicit format override
            result_writer = piru::io::make_result_writer(output_path, output_format);
        } else {
            // Auto-detect from extension
            result_writer = piru::io::make_result_writer(output_path);
        }
        if (!result_writer) {
            LOG_ERROR("map: failed to create output writer for '" + output_path + "'");
            return 1;
        }
        LOG_INFO("Writing results to: " + output_path);
    }

    PIRU_PROFILE_START(profile, "map");

    // =========================================================================
    // INDEX LOADING: Pre-built index OR in-memory indexing
    // =========================================================================
    //
    // Two mutually exclusive workflows:
    //   1. --index: Load pre-serialized index from disk (faster startup)
    //   2. --graph: Run full indexing pipeline in-memory (enables path-walk linearization)
    //
    // Note: Linearization coordinates are only available from in-memory indexing
    //       because path-walk coordinates are not yet serialized to disk.

    std::unique_ptr<piru::index::GraphStore> graph_store;
    std::unique_ptr<piru::index::SeedStore> seed_store;
    std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
    std::vector<std::size_t> path_lengths;  // For anchor bounds checking
    std::string fuzzy_quantizer_name;
    std::string pore_model_name;  // Track pore model for event pipeline parameter selection
    float fuzzy_fine_min{-2.0f};
    float fuzzy_fine_max{2.0f};
    float fuzzy_fine_range{0.4f};
    std::uint32_t fuzzy_n_bins{0};

    PIRU_PROFILE_START(profile, "index");

    if (has_index) {
        // ---------------------------------------------------------------------
        // WORKFLOW 1: Load pre-built index from disk
        // ---------------------------------------------------------------------
        const std::string index_path = parsed.values.at("index");
        LOG_INFO("loading index from " + index_path);

        // Load simple index (.pirx format)
        if (!piru::io::index::is_simple_index(index_path)) {
            LOG_ERROR("map: unsupported index format '" + index_path + "' (expected .pirx file)");
            return 1;
        }

        auto loaded = piru::io::index::load_simple_index(index_path);

        LOG_INFO("simple index loaded: " + std::to_string(loaded.graph->nodeCount()) + " nodes");
        LOG_INFO("index metadata: model=" + loaded.metadata.model_name +
                 ", fuzzy=" + loaded.metadata.fuzzy_quantizer +
                 ", seeds=" + loaded.seeds->extractor_name());

        // Extract path lengths from loaded graph BEFORE moving (AdjListGraphStore has path access)
        const auto& aln_graph = loaded.graph->graph();
        const auto& paths = aln_graph.paths();
        path_lengths.resize(paths.size());
        for (std::size_t i = 0; i < paths.size(); ++i) {
            path_lengths[i] = paths[i].length;
        }

        graph_store = std::move(loaded.graph);
        seed_store = std::move(loaded.seeds);
        linearization_coords = std::move(loaded.linearization_coords);
        fuzzy_quantizer_name = loaded.metadata.fuzzy_quantizer;
        pore_model_name = loaded.metadata.model_name;

    } else {
        // ---------------------------------------------------------------------
        // WORKFLOW 2: Run in-memory indexing
        // ---------------------------------------------------------------------
        const std::string graph_path = parsed.values.at("graph");
        const std::string model_arg = parsed.values.at("model");

        LOG_INFO("Running in-memory indexing: graph=" + graph_path);

        // Step 1: Load pore model (built-in or file)
        auto model = load_model_or_file(model_arg);
        if (!model) {
            return 1;
        }

        // Step 2: Load graph from file
        auto loader = piru::io::make_graph_loader(graph_path);
        if (!loader) {
            LOG_ERROR("map: unsupported graph format for '" + graph_path + "'");
            return 1;
        }

        piru::io::ImportedGraph imported;
        imported.flavor = piru::io::ImportedGraphFlavor::kVg;

        if (!loader->load(imported)) {
            LOG_ERROR("map: failed to read graph file '" + graph_path + "'");
            return 1;
        }

        LOG_INFO("loaded graph: " + std::to_string(imported.nodes.size()) + " nodes, " +
                 std::to_string(imported.edges.size()) + " edges, " +
                 std::to_string(imported.paths.size()) + " paths (" +
                 loader->get_format_name() + ")");

        // Step 3: Configure indexing pipeline
        piru::index::IndexPipelineConfig index_config;

        // Apply CLI overrides for indexer backend
        if (parsed.values.count("indexer-backend")) {
            index_config.indexer_backend = parsed.values.at("indexer-backend");
            LOG_DEBUG("Using indexer backend=" + index_config.indexer_backend);
        }

        // Apply CLI overrides for seed extraction (affects in-memory indexing)
        if (parsed.values.count("seed-k")) {
            index_config.seed_k = std::stoull(parsed.values.at("seed-k"));
            LOG_DEBUG("Using seed k=" + std::to_string(index_config.seed_k) + " for indexing");
        }
        if (parsed.values.count("seed-stride")) {
            index_config.seed_stride = std::stoull(parsed.values.at("seed-stride"));
            LOG_DEBUG("Using seed stride=" + std::to_string(index_config.seed_stride) + " for indexing");
        }
        if (parsed.values.count("seed-backend")) {
            index_config.seed_backend = parsed.values.at("seed-backend");
            LOG_DEBUG("Using seed backend=" + index_config.seed_backend + " for indexing");
        }
        if (parsed.values.count("seed-w")) {
            index_config.seed_window = std::stoull(parsed.values.at("seed-w"));
            LOG_DEBUG("Using seed window=" + std::to_string(index_config.seed_window) + " for indexing");
        }
        if (parsed.values.count("seed-filter")) {
            index_config.seed_filter = std::stod(parsed.values.at("seed-filter"));
            LOG_DEBUG("Using seed filter=" + std::to_string(index_config.seed_filter) + " for indexing");
        }
        if (parsed.values.count("seed-subsample")) {
            index_config.seed_subsample = std::stod(parsed.values.at("seed-subsample"));
            LOG_DEBUG("Using seed subsample=" + std::to_string(index_config.seed_subsample) + " for indexing");
        }
        if (parsed.values.count("seed-mode")) {
            index_config.seed_mode = parsed.values.at("seed-mode");
            LOG_DEBUG("Using seed mode=" + index_config.seed_mode + " for indexing");
        }
        if (parsed.values.count("fuzzy-backend")) {
            index_config.fuzzy_quantizer = parsed.values.at("fuzzy-backend");
            LOG_DEBUG("Using fuzzy quantizer=" + index_config.fuzzy_quantizer + " for indexing");
        }
        if (parsed.values.count("fuzzy-fine-min")) {
            index_config.fuzzy_fine_min = std::stof(parsed.values.at("fuzzy-fine-min"));
            LOG_DEBUG("Using fuzzy fine_min=" + std::to_string(index_config.fuzzy_fine_min));
        }
        if (parsed.values.count("fuzzy-fine-max")) {
            index_config.fuzzy_fine_max = std::stof(parsed.values.at("fuzzy-fine-max"));
            LOG_DEBUG("Using fuzzy fine_max=" + std::to_string(index_config.fuzzy_fine_max));
        }
        if (parsed.values.count("fuzzy-fine-range")) {
            index_config.fuzzy_fine_range = std::stof(parsed.values.at("fuzzy-fine-range"));
            LOG_DEBUG("Using fuzzy fine_range=" + std::to_string(index_config.fuzzy_fine_range));
        }
        if (parsed.values.count("fuzzy-bins")) {
            index_config.fuzzy_n_bins = static_cast<std::uint32_t>(std::stoul(parsed.values.at("fuzzy-bins")));
            LOG_DEBUG("Using fuzzy n_bins=" + std::to_string(index_config.fuzzy_n_bins));
        }

        // Pass executor for parallel indexing
        index_config.executor = executor.get();

        // Step 4: Run full indexing pipeline
        // (simple expand → path-walk index)
        auto index_result = piru::index::run_index_pipeline(imported, *model, index_config);

        graph_store = std::move(index_result.graph_store);
        seed_store = std::move(index_result.seed_store);
        linearization_coords = std::move(index_result.linearization_coords);
        path_lengths = std::move(index_result.path_lengths);
        fuzzy_quantizer_name = index_result.fuzzy_quantizer;
        pore_model_name = model_arg;

        // Preserve fuzzy fine params for query-side processing
        fuzzy_fine_min = index_config.fuzzy_fine_min;
        fuzzy_fine_max = index_config.fuzzy_fine_max;
        fuzzy_fine_range = index_config.fuzzy_fine_range;
        fuzzy_n_bins = index_config.fuzzy_n_bins;

        LOG_INFO("In-memory indexing complete");
    }

    PIRU_PROFILE_STOP(profile, "index");

    // =========================================================================
    // READ FILE DISCOVERY
    // =========================================================================

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

    // =========================================================================
    // MAPPER CONFIGURATION
    // =========================================================================

    piru::mapping::BatchMapperConfig map_config;
    map_config.num_threads = num_threads;
    map_config.seed_store = seed_store.get();
    map_config.graph_store = graph_store.get();
    map_config.linearization_coords = &linearization_coords;  // For DP chaining
    map_config.path_lengths = &path_lengths;  // For anchor bounds checking

    // Configure fuzzy quantizer from index metadata
    if (!fuzzy_quantizer_name.empty()) {
        map_config.fuzzy_config.backend = fuzzy_quantizer_name;
    }
    map_config.fuzzy_config.pore_model = pore_model_name;  // For chemistry-specific defaults
    map_config.fuzzy_config.fine_min = fuzzy_fine_min;
    map_config.fuzzy_config.fine_max = fuzzy_fine_max;
    map_config.fuzzy_config.fine_range = fuzzy_fine_range;
    map_config.fuzzy_config.n_bins = fuzzy_n_bins;
    LOG_DEBUG("Using fuzzy quantizer: " + map_config.fuzzy_config.backend +
              " (fine_min=" + std::to_string(fuzzy_fine_min) +
              ", fine_max=" + std::to_string(fuzzy_fine_max) +
              ", fine_range=" + std::to_string(fuzzy_fine_range) +
              ", n_bins=" + std::to_string(fuzzy_n_bins) + ")");

    // Configure seed extractor from seed store parameters
    const auto* hash_seed_store = dynamic_cast<piru::index::HashSeedStore*>(seed_store.get());
    if (!hash_seed_store) {
        LOG_ERROR("map: unsupported seed store backend");
        return 1;
    }

    auto index_seed_cfg = seed_config_from_store(*hash_seed_store);
    if (index_seed_cfg.backend.empty()) {
        LOG_ERROR("map: seed store missing extractor name");
        return 1;
    }

    LOG_DEBUG("Using seed extractor settings: backend=" + index_seed_cfg.backend +
              ", k=" + std::to_string(index_seed_cfg.k) +
              ", stride=" + std::to_string(index_seed_cfg.stride) +
              ", window=" + std::to_string(index_seed_cfg.window) +
              ", qbits=" + std::to_string(index_seed_cfg.qbits));
    map_config.seed_config = index_seed_cfg;

    // Configure clusterer from seed store statistics
    map_config.clusterer_config.max_hash_frequency = hash_seed_store->max_hash_frequency();
    LOG_DEBUG("Using clustering config from index: max_hash_frequency=" +
              std::to_string(map_config.clusterer_config.max_hash_frequency));

    // Configure clusterer backend (default: dp-chain for path-walk pipeline)
    const std::string clusterer = parsed.values.count("clusterer")
                                   ? parsed.values.at("clusterer")
                                   : "dp-chain";
    map_config.clusterer_config.backend = clusterer;

    // DP chain parameter overrides (only relevant when clusterer=dp-chain)
    if (parsed.values.count("chain-max-dist")) {
        map_config.clusterer_config.dp_max_dist = std::stoull(parsed.values.at("chain-max-dist"));
    }
    if (parsed.values.count("chain-max-diag-dev")) {
        map_config.clusterer_config.dp_max_diag_dev = std::stoull(parsed.values.at("chain-max-diag-dev"));
    }
    if (parsed.values.count("chain-gap-penalty")) {
        map_config.clusterer_config.dp_gap_penalty = std::stod(parsed.values.at("chain-gap-penalty"));
    }
    if (parsed.values.count("chain-diag-penalty")) {
        map_config.clusterer_config.dp_diag_penalty = std::stod(parsed.values.at("chain-diag-penalty"));
    }
    if (parsed.values.count("chain-overlap-penalty")) {
        map_config.clusterer_config.dp_overlap_penalty = std::stod(parsed.values.at("chain-overlap-penalty"));
    }
    if (parsed.values.count("chain-anchor-weight")) {
        map_config.clusterer_config.dp_anchor_weight = std::stod(parsed.values.at("chain-anchor-weight"));
    }
    if (parsed.values.count("chain-min-score")) {
        map_config.clusterer_config.dp_min_chain_score = std::stoull(parsed.values.at("chain-min-score"));
    }
    if (parsed.values.count("chain-max-chains")) {
        map_config.clusterer_config.dp_max_chains = std::stoull(parsed.values.at("chain-max-chains"));
    }
    if (parsed.values.count("chain-max-skip")) {
        map_config.clusterer_config.dp_max_skip = std::stoull(parsed.values.at("chain-max-skip"));
    }
    if (parsed.values.count("chain-merge")) {
        const std::string val = parsed.values.at("chain-merge");
        map_config.clusterer_config.dp_merge_chains = (val == "true" || val == "1" || val == "yes");
    }

    // Configure event pipeline (unified event detection + normalization)
    map_config.event_pipeline_config.pore_model = pore_model_name;
    if (parsed.values.count("event-pipeline")) {
        map_config.event_pipeline_config.backend = parsed.values.at("event-pipeline");
    }
    // Event detection parameter overrides (take precedence over backend defaults)
    if (parsed.values.count("event-w1")) {
        map_config.event_pipeline_config.override_window_length1 = std::stoi(parsed.values.at("event-w1"));
    }
    if (parsed.values.count("event-w2")) {
        map_config.event_pipeline_config.override_window_length2 = std::stoi(parsed.values.at("event-w2"));
    }
    if (parsed.values.count("event-t1")) {
        map_config.event_pipeline_config.override_threshold1 = std::stof(parsed.values.at("event-t1"));
    }
    if (parsed.values.count("event-t2")) {
        map_config.event_pipeline_config.override_threshold2 = std::stof(parsed.values.at("event-t2"));
    }
    if (parsed.values.count("event-peak")) {
        map_config.event_pipeline_config.override_peak_height = std::stof(parsed.values.at("event-peak"));
    }

    // Map-time frequency cap: recompute threshold from percentile
    if (parsed.values.count("seed-freq-cap")) {
        const double percentile = std::stod(parsed.values.at("seed-freq-cap"));
        const_cast<piru::index::HashSeedStore*>(hash_seed_store)->recompute_threshold_from_percentile(percentile);
        LOG_INFO("Map-time seed freq cap at " + std::to_string(percentile * 100) + "th percentile → threshold=" +
                 std::to_string(hash_seed_store->frequency_threshold()));
    }

    // Log the seed frequency threshold being used
    LOG_INFO("Seed frequency threshold: " + std::to_string(hash_seed_store->frequency_threshold()) +
             " (seeds with freq > " + std::to_string(hash_seed_store->frequency_threshold()) + " will be skipped)");

    // Dump seed store if requested (do this before mapping starts)
    if (parsed.values.count("dump-seed-store")) {
        const std::string dump_path = parsed.values.at("dump-seed-store");
        std::ofstream out(dump_path);
        if (!out.is_open()) {
            LOG_ERROR("Failed to open seed store dump file: " + dump_path);
            return 1;
        }
        out << "# Seed store dump: hash -> frequency\n";
        out << "# freq_threshold=" << hash_seed_store->frequency_threshold() << "\n";
        out << "hash\tfrequency\n";
        for (const auto& [hash, hits] : hash_seed_store->data()) {
            out << std::hex << hash << std::dec << "\t" << hits.size() << "\n";
        }
        out.close();
        LOG_INFO("Dumped seed store (" + std::to_string(hash_seed_store->size()) + " hashes) to: " + dump_path);
    }

    // Configure result writer if specified
    if (result_writer) {
        map_config.result_writer = result_writer.get();
    }

    // Configure debug dump directories
    if (parsed.values.count("dump-anchors")) {
        map_config.dump_anchors_dir = parsed.values.at("dump-anchors");
        std::filesystem::create_directories(map_config.dump_anchors_dir);
        LOG_INFO("Dumping anchors to: " + map_config.dump_anchors_dir);
    }
    if (parsed.values.count("dump-chains")) {
        map_config.dump_chains_dir = parsed.values.at("dump-chains");
        std::filesystem::create_directories(map_config.dump_chains_dir);
        LOG_INFO("Dumping chains to: " + map_config.dump_chains_dir);
    }
    if (parsed.values.count("dump-hit-stats")) {
        map_config.dump_hit_stats_dir = parsed.values.at("dump-hit-stats");
        std::filesystem::create_directories(map_config.dump_hit_stats_dir);
        LOG_INFO("Dumping hit stats to: " + map_config.dump_hit_stats_dir);
    }
    if (parsed.values.count("dump-path-chains")) {
        map_config.dump_path_chains_dir = parsed.values.at("dump-path-chains");
        std::filesystem::create_directories(map_config.dump_path_chains_dir);
        LOG_INFO("Dumping path chains to: " + map_config.dump_path_chains_dir);
    }
    if (parsed.values.count("dump-read-seeds")) {
        map_config.dump_read_seeds_dir = parsed.values.at("dump-read-seeds");
        std::filesystem::create_directories(map_config.dump_read_seeds_dir);
        LOG_INFO("Dumping read seeds to: " + map_config.dump_read_seeds_dir);
    }
    if (parsed.values.count("dump-anchor-detail")) {
        map_config.dump_anchor_detail_dir = parsed.values.at("dump-anchor-detail");
        std::filesystem::create_directories(map_config.dump_anchor_detail_dir);
        LOG_INFO("Dumping anchor detail to: " + map_config.dump_anchor_detail_dir);
    }
    if (parsed.values.count("dump-anchor-reads")) {
        std::string reads_str = parsed.values.at("dump-anchor-reads");
        std::istringstream iss(reads_str);
        std::string token;
        while (std::getline(iss, token, ',')) {
            if (!token.empty()) {
                map_config.dump_anchor_reads.insert(token);
            }
        }
        LOG_INFO("Dumping anchor detail for " + std::to_string(map_config.dump_anchor_reads.size()) + " read filters");
    }
    if (parsed.values.count("no-anchor-merge")) {
        map_config.enable_anchor_merge = false;
        LOG_INFO("Anchor merging disabled");
    }

    // Configure result formatter
    if (parsed.values.count("min-secondary-ratio")) {
        map_config.formatter_config.min_secondary_ratio =
            std::stod(parsed.values.at("min-secondary-ratio"));
    }

    // =========================================================================
    // READ PROCESSING
    // =========================================================================

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
             ", batches=" + std::to_string(total_batches) +
             ", reads=" + std::to_string(total_reads));

    PIRU_PROFILE_STOP(profile, "map");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
