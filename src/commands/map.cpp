#include "commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/index_pipeline.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/index/serialization.hpp"
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
        {'\0', "", false, "\nIn-Memory Indexing Options (with --graph):"},
        {'\0', "linearizer", true, "Linearizer backend: superbubble (default) or path-walk"},
        {'\0', "graph-type", true, "Graph type: dbg (default) or vg"},
        {'\0', "graph-k", true, "DBG k-mer size (default: auto-detect from overlap)"},
        {'\0', "", false, "\nMapping Options:"},
        {'\0', "max-seed-freq", true, "Maximum seed frequency for lookup (default: use index threshold)"},
        {'\0', "clusterer", true, "Clusterer backend: fse (default), probe, dp-chain"},
        {'\0', "align", false, "Enable signal-level alignment for chain evaluation"},
        {'\0', "align-backend", true, "Alignment backend: path-guided (default), radius, auto"},
        {'\0', "", false, "\nSignal Processing Options (only with --graph):"},
        {'\0', "event-pipeline", true, "Event pipeline backend: scrappie (default), rawhash, passthrough"},
        {'\0', "fuzzy-backend", true, "Fuzzy quantizer backend (default: rh2)"},
        {'\0', "fuzzy-fine-min", true, "Fuzzy quantizer fine region min (default: -2.0)"},
        {'\0', "fuzzy-fine-max", true, "Fuzzy quantizer fine region max (default: 2.0)"},
        {'\0', "fuzzy-fine-range", true, "Fuzzy quantizer fine bin width (default: 0.4)"},
        {'\0', "seed-k", true, "Seed extractor k-mer size (default: 6)"},
        {'\0', "seed-stride", true, "Seed extractor stride (default: 1)"},
        {'\0', "", false, "\nDebug Options:"},
        {'\0', "dump-anchors", true, "Dump anchors to directory (one file per read)"},
        {'\0', "dump-chains", true, "Dump chains to directory (one file per read)"},
        {'\0', "", false, "\nOutput Options:"},
        {'o', "output", true, "Output file path (format auto-detected from extension: .paf, .gaf, .gam, .json)"},
        {'\0', "output-format", true, "Override output format (paf, gaf, gam, json)"},
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

    // Extract common options
    const bool profile = parsed.values.count("profile") > 0;
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
    std::unique_ptr<piru::index::SignalStore> signal_store;
    std::unique_ptr<piru::index::SeedStore> seed_store;
    std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
    std::string fuzzy_quantizer_name;
    std::string pore_model_name;  // Track pore model for event pipeline parameter selection
    float fuzzy_fine_min{-2.0f};
    float fuzzy_fine_max{2.0f};
    float fuzzy_fine_range{0.4f};

    if (has_index) {
        // ---------------------------------------------------------------------
        // WORKFLOW 1: Load pre-built index from disk
        // ---------------------------------------------------------------------
        const std::string index_path = parsed.values.at("index");
        LOG_INFO("loading index from " + index_path);

        auto loaded_index = piru::io::index::load_index(index_path);
        if (!loaded_index.graph) {
            LOG_ERROR("map: failed to load index from '" + index_path + "'");
            return 1;
        }
        if (!loaded_index.seeds) {
            LOG_ERROR("map: index is missing seeds; cannot perform lookups");
            return 1;
        }

        LOG_INFO("index loaded: " + std::to_string(loaded_index.graph->nodeCount()) + " nodes");
        LOG_INFO("index metadata: fuzzy=" + loaded_index.metadata.fuzzy_quantizer +
                 ", align=" + loaded_index.metadata.align_quantizer +
                 ", seeds=" + loaded_index.seeds->extractor_name());

        graph_store = std::move(loaded_index.graph);
        signal_store = std::move(loaded_index.signals);
        seed_store = std::move(loaded_index.seeds);
        fuzzy_quantizer_name = loaded_index.metadata.fuzzy_quantizer;
        pore_model_name = loaded_index.metadata.model_name;

        // Note: Linearization coords not serialized yet (Phase 7 TODO)
        // Superbubble uses graph.node.chain_id/linear_position instead

        // Warn if user tries to override seed/fuzzy params with pre-built index
        if (parsed.values.count("seed-k") || parsed.values.count("seed-stride") ||
            parsed.values.count("fuzzy-backend") || parsed.values.count("fuzzy-fine-min") ||
            parsed.values.count("fuzzy-fine-max") || parsed.values.count("fuzzy-fine-range")) {
            LOG_WARN("--seed-* and --fuzzy-* flags are ignored with --index (use --graph for experimentation)");
        }

    } else {
        // ---------------------------------------------------------------------
        // WORKFLOW 2: Run in-memory indexing
        // ---------------------------------------------------------------------
        const std::string graph_path = parsed.values.at("graph");
        const std::string model_arg = parsed.values.at("model");
        const std::string linearizer = parsed.values.count("linearizer")
                                           ? parsed.values.at("linearizer")
                                           : "superbubble";
        const std::string graph_type = parsed.values.count("graph-type")
                                           ? parsed.values.at("graph-type")
                                           : "dbg";

        LOG_INFO("running in-memory indexing: graph=" + graph_path +
                 ", linearizer=" + linearizer);

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
        if (graph_type == "vg") {
            imported.flavor = piru::io::ImportedGraphFlavor::kVg;
        } else if (graph_type == "dbg") {
            imported.flavor = piru::io::ImportedGraphFlavor::kDbg;
        } else {
            imported.flavor = piru::io::ImportedGraphFlavor::kUnknown;
        }

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
        index_config.graph_flavor = graph_type;
        index_config.linearizer = linearizer;

        // Auto-detect graph_k from edge overlap (DBG only)
        if (!parsed.values.count("graph-k") && graph_type == "dbg") {
            if (!imported.edges.empty() && imported.edges.front().overlap_bases.has_value()) {
                const auto overlap = *imported.edges.front().overlap_bases;
                index_config.graph_k = overlap + 1;
            }
        } else if (parsed.values.count("graph-k")) {
            index_config.graph_k = std::stoul(parsed.values.at("graph-k"));
        }

        // Apply CLI overrides for seed extraction (affects in-memory indexing)
        if (parsed.values.count("seed-k")) {
            index_config.seed_k = std::stoull(parsed.values.at("seed-k"));
            LOG_INFO("Using seed k=" + std::to_string(index_config.seed_k) + " for indexing");
        }
        if (parsed.values.count("seed-stride")) {
            index_config.seed_stride = std::stoull(parsed.values.at("seed-stride"));
            LOG_INFO("Using seed stride=" + std::to_string(index_config.seed_stride) + " for indexing");
        }
        if (parsed.values.count("fuzzy-backend")) {
            index_config.fuzzy_quantizer = parsed.values.at("fuzzy-backend");
            LOG_INFO("Using fuzzy quantizer=" + index_config.fuzzy_quantizer + " for indexing");
        }
        if (parsed.values.count("fuzzy-fine-min")) {
            index_config.fuzzy_fine_min = std::stof(parsed.values.at("fuzzy-fine-min"));
            LOG_INFO("Using fuzzy fine_min=" + std::to_string(index_config.fuzzy_fine_min));
        }
        if (parsed.values.count("fuzzy-fine-max")) {
            index_config.fuzzy_fine_max = std::stof(parsed.values.at("fuzzy-fine-max"));
            LOG_INFO("Using fuzzy fine_max=" + std::to_string(index_config.fuzzy_fine_max));
        }
        if (parsed.values.count("fuzzy-fine-range")) {
            index_config.fuzzy_fine_range = std::stof(parsed.values.at("fuzzy-fine-range"));
            LOG_INFO("Using fuzzy fine_range=" + std::to_string(index_config.fuzzy_fine_range));
        }

        // Step 4: Run full indexing pipeline
        // (transform → linearize → squigglize → quantize → seed extraction)
        auto index_result = piru::index::run_index_pipeline(imported, *model, index_config);

        graph_store = std::move(index_result.graph_store);
        signal_store = std::move(index_result.signal_store);
        seed_store = std::move(index_result.seed_store);
        linearization_coords = std::move(index_result.linearization_coords);
        fuzzy_quantizer_name = index_result.fuzzy_quantizer;
        pore_model_name = model_arg;

        // Preserve fuzzy fine params for query-side processing
        fuzzy_fine_min = index_config.fuzzy_fine_min;
        fuzzy_fine_max = index_config.fuzzy_fine_max;
        fuzzy_fine_range = index_config.fuzzy_fine_range;

        LOG_INFO("in-memory indexing complete");
    }

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

    // Configure fuzzy quantizer from index metadata
    if (!fuzzy_quantizer_name.empty()) {
        map_config.fuzzy_config.backend = fuzzy_quantizer_name;
    }
    map_config.fuzzy_config.fine_min = fuzzy_fine_min;
    map_config.fuzzy_config.fine_max = fuzzy_fine_max;
    map_config.fuzzy_config.fine_range = fuzzy_fine_range;
    LOG_INFO("Using fuzzy quantizer: " + map_config.fuzzy_config.backend +
             " (fine_min=" + std::to_string(fuzzy_fine_min) +
             ", fine_max=" + std::to_string(fuzzy_fine_max) +
             ", fine_range=" + std::to_string(fuzzy_fine_range) + ")");

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

    LOG_INFO("Using seed extractor settings: backend=" + index_seed_cfg.backend +
             ", k=" + std::to_string(index_seed_cfg.k) +
             ", stride=" + std::to_string(index_seed_cfg.stride) +
             ", qbits=" + std::to_string(index_seed_cfg.qbits));
    map_config.seed_config = index_seed_cfg;

    // Configure clusterer from seed store statistics
    map_config.clusterer_config.max_hash_frequency = hash_seed_store->max_hash_frequency();
    LOG_INFO("Using clustering config from index: max_hash_frequency=" +
             std::to_string(map_config.clusterer_config.max_hash_frequency));

    // Configure clusterer backend (default: fse for backward compatibility)
    const std::string clusterer = parsed.values.count("clusterer")
                                   ? parsed.values.at("clusterer")
                                   : "fse";
    map_config.clusterer_config.backend = clusterer;

    // Configure event pipeline (unified event detection + normalization)
    map_config.event_pipeline_config.pore_model = pore_model_name;
    if (parsed.values.count("event-pipeline")) {
        map_config.event_pipeline_config.backend = parsed.values.at("event-pipeline");
    }
    // Default backend is "scrappie" (set in EventPipelineConfig)

    // Override seed frequency threshold if specified
    if (parsed.values.count("max-seed-freq")) {
        const std::size_t override_threshold = std::stoull(parsed.values.at("max-seed-freq"));
        const_cast<piru::index::HashSeedStore*>(hash_seed_store)->set_frequency_threshold(override_threshold);
        LOG_INFO("Overriding seed frequency threshold to " + std::to_string(override_threshold));
    }

    // Configure result writer if specified
    if (result_writer) {
        map_config.result_writer = result_writer.get();
    }

    // Configure signal store for alignment (if available)
    map_config.signal_store = signal_store.get();

    // Configure alignment if enabled
    if (parsed.values.count("align")) {
        map_config.enable_alignment = true;

        // Parse alignment backend
        const std::string align_backend = parsed.values.count("align-backend")
                                              ? parsed.values.at("align-backend")
                                              : "path-guided";
        if (align_backend == "radius") {
            map_config.align_config.backend = piru::alignment::AlignerBackend::kRadius;
        } else if (align_backend == "auto") {
            map_config.align_config.backend = piru::alignment::AlignerBackend::kAuto;
        } else {
            map_config.align_config.backend = piru::alignment::AlignerBackend::kPathGuided;
        }

        LOG_INFO("Signal-level alignment enabled: backend=" + align_backend);
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

    // =========================================================================
    // READ PROCESSING
    // =========================================================================

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
        LOG_INFO("map: processing '" + f.string() + "' (format=" + provider->get_format_name() +
                 ").");

        piru::mapping::BatchMapper mapper(*provider, map_config, std::cout);
        const auto stats = mapper.process_all();
        total_reads += stats.reads_processed;
        total_batches += stats.batches;
    }

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
