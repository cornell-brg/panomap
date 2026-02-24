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
#include "io/index/simple_serialization.hpp"
#include "io/models/model_factory.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

#ifdef PIRU_DUMP_GRAPHS
#include "io/graphs/gfa_exporter.hpp"
#endif

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
        {'m', "model", true, "Pore model (builtin name: r9.4/r10.4 or model file path)"},
        {'o', "output", true, "Output index file (default: <graph-file>.pirx)"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
        {'\0', "", false, "\nSeed Generation Options:"},
        {'\0', "seed-backend", true, "Seed extractor backend: kmer (default), minimizer"},
        {'\0', "seed-k", true, "Seed k-mer size (default: 6)"},
        {'\0', "seed-w", true, "Minimizer window size (default: 5, only with --seed-backend minimizer)"},
        {'\0', "seed-stride", true, "Seed stride (default: 1)"},
        {'\0', "seed-filter", true, "Keep least frequent seed fraction (default: 1.0)"},
        {'\0', "seed-mode", true, "Seeding mode: node, path (default)"},
        {'\0', "", false, "\nIndexer Options:"},
        {'\0', "indexer-backend", true, "Indexer backend: node-first, path-walk (default)"},
        {'\0', "", false, "\nDebug Options:"},
        {'\0', "dump-norm-stats", true, "Dump per-path normalization stats to TSV file (path-walk only)"},
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

    // Parse thread count and create executor
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
    imported.flavor = piru::io::ImportedGraphFlavor::kVg;

    if (!loader->load(imported)) {
        LOG_ERROR("index: failed to read graph file '" + graph_path + "'");
        return 1;
    }

    LOG_INFO("loaded graph: " + std::to_string(imported.nodes.size()) + " nodes, " +
             std::to_string(imported.edges.size()) + " edges, " +
             std::to_string(imported.paths.size()) + " paths (" + loader->get_format_name() + ")");

#ifdef PIRU_DUMP_GRAPHS
    piru::GfaExporter::dumpImportedGraph(imported, "imported_graph.gfa");
#endif

    // -------------------------------------------------------------------------
    // Parameter parsing and validation
    // -------------------------------------------------------------------------

    const std::size_t pore_k = model->k();

    // Use defaults from IndexPipelineConfig (single source of truth)
    piru::index::IndexPipelineConfig defaults;
    const std::string seed_backend =
        parsed.values.count("seed-backend") ? parsed.values.at("seed-backend") : defaults.seed_backend;
    const std::size_t seed_k =
        parsed.values.count("seed-k") ? std::stoul(parsed.values.at("seed-k")) : defaults.seed_k;
    const std::size_t seed_w =
        parsed.values.count("seed-w") ? std::stoul(parsed.values.at("seed-w")) : defaults.seed_window;
    const std::size_t seed_stride =
        parsed.values.count("seed-stride") ? std::stoul(parsed.values.at("seed-stride")) : defaults.seed_stride;
    const double seed_filter = parsed.values.count("seed-filter")
                                   ? std::stod(parsed.values.at("seed-filter"))
                                   : defaults.seed_filter;
    const std::string seed_mode = parsed.values.count("seed-mode")
                                      ? parsed.values.at("seed-mode")
                                      : defaults.seed_mode;

    // Default output in current directory (stem without extension)
    std::string output_base = std::filesystem::path(graph_path).stem().string();
    if (parsed.values.count("output")) {
        output_base = parsed.values.at("output");
    }

    LOG_INFO("input: " + graph_path);
    LOG_INFO("model: " + model->name() + " (k=" + std::to_string(pore_k) + ")");
    LOG_INFO("seeds: backend=" + seed_backend + ", k=" + std::to_string(seed_k) +
             ", w=" + std::to_string(seed_w) + ", stride=" + std::to_string(seed_stride) +
             ", filter=" + std::to_string(seed_filter) + ", mode=" + seed_mode);
    LOG_INFO("output: " + output_base);

    // -------------------------------------------------------------------------
    // Run Indexing Pipeline (shared with map --graph)
    // -------------------------------------------------------------------------

    piru::index::IndexPipelineConfig index_config;
    index_config.seed_backend = seed_backend;
    index_config.seed_k = seed_k;
    index_config.seed_window = seed_w;
    index_config.seed_stride = seed_stride;
    index_config.seed_filter = seed_filter;
    index_config.seed_mode = seed_mode;
    index_config.fuzzy_quantizer = "rh2";
    if (parsed.values.count("indexer-backend")) {
        index_config.indexer_backend = parsed.values.at("indexer-backend");
    }
    if (parsed.values.count("dump-norm-stats")) {
        index_config.dump_norm_stats_path = parsed.values.at("dump-norm-stats");
    }

    // Pass executor for parallel indexing
    index_config.executor = executor.get();

    auto result = piru::index::run_index_pipeline(imported, *model, index_config);

    // -------------------------------------------------------------------------
    // Serialize Index
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "serialize");

    std::string output_path = output_base;
    if (!output_path.ends_with(".pirx")) {
        output_path += ".pirx";
    }

    piru::io::index::SimpleIndexMetadata metadata;
    metadata.model_name = result.model_name;
    metadata.pore_k = result.pore_k;
    metadata.fuzzy_quantizer = result.fuzzy_quantizer;
    metadata.graph_flavor = "vg";

    piru::io::index::save_simple_index(
        output_path,
        *result.graph_store,
        *result.seed_store,
        result.linearization_coords,
        metadata);

    LOG_INFO("Index written to " + output_path);

    PIRU_PROFILE_STOP(profile, "serialize");

    PIRU_PROFILE_STOP(profile, "index");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
