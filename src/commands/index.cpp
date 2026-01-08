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
        {'\0', "index-backend", true, "Index pipeline: classic (default) or simple"},
        {'g', "graph", true, "Graph type: vg (default) or dbg"},
        {'v', "verbose", false, "Enable verbose logging (DEBUG level)"},
        {'k', "graph-k", true, "DBG k-mer size (default: auto-detect from overlap)"},
        {'m', "model", true, "Pore model (builtin name: r9.4/r10.4 or model file path)"},
        {'o', "output", true, "Output index directory (default: <graph-file>.piru)"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
        {'\0', "", false, "\nSeed Generation Options:"},
        {'\0', "seed-k", true, "Seed k-mer size (default: 6)"},
        {'\0', "seed-stride", true, "Seed stride (default: 1)"},
        {'\0', "seed-filter", true, "Keep least frequent seed fraction (default: 0.9)"},
        {'\0', "seed-mode", true, "Seeding mode: node, path (default)"},
        {'\0', "", false, "\nAlignment Quantization Options:"},
        {'\0', "aq-backend", true, "Backend: int16 (default), int8, passthrough"},
        {'\0', "aq-scale", true, "Manual scale override (expert)"},
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
    const std::string model_arg = parsed.values.count("model") ? parsed.values.at("model") : "r10.4";
    const std::string graph_type = parsed.values.count("graph") ? parsed.values.at("graph") : "vg";

    if (parsed.positionals.empty()) {
        LOG_ERROR("index: missing required <graph-file>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    const std::string graph_path = parsed.positionals.front();
    const std::string graph_flavor = graph_type;

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
    if (graph_flavor == "vg") {
        imported.flavor = piru::io::ImportedGraphFlavor::kVg;
    } else if (graph_flavor == "dbg") {
        imported.flavor = piru::io::ImportedGraphFlavor::kDbg;
    } else {
        imported.flavor = piru::io::ImportedGraphFlavor::kUnknown;
    }

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
    std::size_t graph_k = 0;

    // Auto-detect graph_k from overlap (DBG case)
    if (!parsed.values.count("graph-k") && graph_flavor == "dbg") {
        if (!imported.edges.empty() && imported.edges.front().overlap_bases.has_value()) {
            const auto overlap = *imported.edges.front().overlap_bases;
            graph_k = overlap + 1;
        }
    } else if (parsed.values.count("graph-k")) {
        graph_k = std::stoul(parsed.values.at("graph-k"));
    }

    // VG graphs don't need graph_k
    if (graph_flavor != "vg") {
        if (graph_k == 0) {
            LOG_ERROR("index: cannot determine graph k-mer size (provide --graph-k)");
            return 1;
        }

        if (graph_k < pore_k) {
            LOG_ERROR("index: graph k=" + std::to_string(graph_k) + " < pore k=" +
                      std::to_string(pore_k) + " (invalid)");
            return 1;
        }
    }

    // Use defaults from IndexPipelineConfig (single source of truth)
    piru::index::IndexPipelineConfig defaults;
    const std::size_t seed_k =
        parsed.values.count("seed-k") ? std::stoul(parsed.values.at("seed-k")) : defaults.seed_k;
    const std::size_t seed_stride =
        parsed.values.count("seed-stride") ? std::stoul(parsed.values.at("seed-stride")) : defaults.seed_stride;
    const double seed_filter = parsed.values.count("seed-filter")
                                   ? std::stod(parsed.values.at("seed-filter"))
                                   : defaults.seed_filter;
    const std::string seed_mode = parsed.values.count("seed-mode")
                                      ? parsed.values.at("seed-mode")
                                      : defaults.seed_mode;

    std::string output_dir = graph_path + ".piru";
    if (parsed.values.count("output")) {
        output_dir = parsed.values.at("output");
    }

    LOG_INFO("input: " + graph_path + " (type=" + graph_flavor + ", graph_k=" +
             std::to_string(graph_k) + ")");
    LOG_INFO("model: " + model->name() + " (k=" + std::to_string(pore_k) + ")");
    LOG_INFO("seeds: k=" + std::to_string(seed_k) + ", stride=" + std::to_string(seed_stride) +
             ", filter=" + std::to_string(seed_filter) + ", mode=" + seed_mode);
    LOG_INFO("output: " + output_dir);

    // -------------------------------------------------------------------------
    // Run Indexing Pipeline (shared with map --graph)
    // -------------------------------------------------------------------------

    piru::index::IndexPipelineConfig index_config;
    index_config.graph_flavor = graph_flavor;
    index_config.linearizer = "superbubble";  // index command always uses superbubble
    if (parsed.values.count("index-backend")) {
        index_config.pipeline_mode = parsed.values.at("index-backend");
        LOG_INFO("Using index backend: " + index_config.pipeline_mode);
    }
    index_config.graph_k = graph_k;
    index_config.seed_k = seed_k;
    index_config.seed_stride = seed_stride;
    index_config.seed_filter = seed_filter;
    index_config.seed_mode = seed_mode;
    index_config.fuzzy_quantizer = "rh2";
    index_config.alignment_quantizer =
        parsed.values.count("aq-backend") ? parsed.values.at("aq-backend") : "int16";
    if (parsed.values.count("aq-scale")) {
        index_config.alignment_scale = std::stod(parsed.values.at("aq-scale"));
    }

    auto result = piru::index::run_index_pipeline(imported, *model, index_config);

    // -------------------------------------------------------------------------
    // Serialize Index
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "serialize");

    if (index_config.pipeline_mode == "simple") {
        // Simple index: single file format (.pir2)
        std::string output_path = output_dir;
        if (!output_path.ends_with(".pir2")) {
            output_path += ".pir2";
        }

        piru::io::index::SimpleIndexMetadata metadata;
        metadata.model_name = result.model_name;
        metadata.pore_k = result.pore_k;
        metadata.fuzzy_quantizer = result.fuzzy_quantizer;
        metadata.graph_flavor = graph_flavor;

        piru::io::index::save_simple_index(
            output_path,
            *result.graph_store,
            *result.seed_store,
            result.linearization_coords,
            metadata);

        LOG_INFO("Simple index written to " + output_path);
    } else {
        // Classic index: directory with multiple files
        if (!result.signal_store) {
            LOG_ERROR("Classic pipeline requires signal store");
            return 1;
        }

        piru::io::index::IndexMetadata metadata;
        metadata.graph_flavor = result.graph_flavor == piru::io::ImportedGraphFlavor::kDbg ? 1 : 2;
        metadata.graph_k = result.graph_k;
        metadata.pore_k = result.pore_k;
        metadata.model_name = result.model_name;
        metadata.fuzzy_quantizer = result.fuzzy_quantizer;
        metadata.align_quantizer = result.alignment_quantizer;
        metadata.source_path = graph_path;

        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directory(output_dir);
        }

        std::filesystem::path dir_path(output_dir);
        std::string basename = dir_path.filename().string();

        piru::io::index::write_graph(dir_path / (basename + ".graph"), *result.graph_store, metadata);
        piru::io::index::write_signals(dir_path / (basename + ".signals"), *result.signal_store,
                                        result.alignment_scale, result.alignment_offset);
        piru::io::index::write_seeds(dir_path / (basename + ".seeds"), *result.seed_store);

        LOG_INFO("Index written to " + output_dir);
    }

    PIRU_PROFILE_STOP(profile, "serialize");

    PIRU_PROFILE_STOP(profile, "index");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
