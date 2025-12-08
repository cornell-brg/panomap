#include "commands/index.hpp"

#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/pseudo_linearize.hpp"
#include "index/seed_builder.hpp"
#include "index/squigglize.hpp"
#include "index/transform_dbg.hpp"
#include "io/graphs/graph.hpp"
#include "io/graphs/graph_loader_factory.hpp"
#include "io/index/serialization.hpp"
#include "io/models/model_factory.hpp"
#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
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
        {'g', "graph", true, "Graph type: dbg (default) or vg"},
        {'k', "graph-k", true, "DBG k-mer size (default: auto-detect from overlap)"},
        {'m', "model", true, "Pore model (builtin name: r9.4/r10.4 or model file path)"},
        {'o', "output", true, "Output index directory (default: <graph-file>.piru)"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
        {'\0', "seed-k", true, "Seed k-mer size (default: 10)"},
        {'\0', "seed-stride", true, "Seed stride (default: 1)"},
        {'\0', "seed-filter", true, "Keep least frequent seed fraction (default: 1.0)"},
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
    const std::string model_arg = parsed.values.count("model") ? parsed.values.at("model") : "r10.4";
    const std::string graph_type = parsed.values.count("graph") ? parsed.values.at("graph") : "dbg";

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

    if (graph_k == 0) {
        LOG_ERROR("index: cannot determine graph k-mer size (provide --graph-k)");
        return 1;
    }

    if (graph_k < pore_k) {
        LOG_ERROR("index: graph k=" + std::to_string(graph_k) + " < pore k=" +
                  std::to_string(pore_k) + " (invalid)");
        return 1;
    }

    const std::size_t seed_k =
        parsed.values.count("seed-k") ? std::stoul(parsed.values.at("seed-k")) : 10;
    const std::size_t seed_stride =
        parsed.values.count("seed-stride") ? std::stoul(parsed.values.at("seed-stride")) : 1;
    const double seed_filter = parsed.values.count("seed-filter")
                                   ? std::stod(parsed.values.at("seed-filter"))
                                   : 1.0;

    std::string output_dir = graph_path + ".piru";
    if (parsed.values.count("output")) {
        output_dir = parsed.values.at("output");
    }

    LOG_INFO("input: " + graph_path + " (type=" + graph_flavor + ", graph_k=" +
             std::to_string(graph_k) + ")");
    LOG_INFO("model: " + model->name() + " (k=" + std::to_string(pore_k) + ")");
    LOG_INFO("seeds: k=" + std::to_string(seed_k) + ", stride=" + std::to_string(seed_stride) +
             ", filter=" + std::to_string(seed_filter));
    LOG_INFO("output: " + output_dir);

    // -------------------------------------------------------------------------
    // Stage 1: Graph Transformation (ImportedGraph -> AlnGraph)
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "transform");
    piru::index::AlnGraph aln_graph;

    if (graph_flavor == "dbg") {
        aln_graph = piru::index::transformDbg(imported, graph_k, pore_k);
    } else if (graph_flavor == "vg") {
        LOG_ERROR("index: VG transformation not yet implemented");
        return 1;
    } else {
        LOG_ERROR("index: unknown graph flavor: " + graph_flavor);
        return 1;
    }

    if (!aln_graph.validate()) {
        LOG_ERROR("index: AlnGraph validation failed after transformation");
        return 1;
    }

    LOG_INFO("[1/4] transformed: " + std::to_string(aln_graph.nodeCount()) + " directional nodes");
    PIRU_PROFILE_STOP(profile, "transform");

    // -------------------------------------------------------------------------
    // Stage 2: Pseudo-Linearization
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "linearize");

    auto scc = piru::index::computeScc(aln_graph);
    auto tips = piru::index::chainTips(aln_graph, scc);
    piru::index::chainCycles(aln_graph, tips);
    auto sb = piru::index::chainSuperbubbles(aln_graph, scc, tips);
    auto chain_ids = piru::index::assignChainIds(sb.uf);
    auto positions = piru::index::assignLinearPositions(aln_graph, chain_ids, scc);

    // Store chain metadata in graph nodes
    for (std::size_t i = 0; i < aln_graph.nodeCount(); ++i) {
        aln_graph.mutableNode(i).chain_id = chain_ids[i];
        aln_graph.mutableNode(i).linear_position = positions[i];
    }

    std::set<std::size_t> unique_chains(chain_ids.begin(), chain_ids.end());
    LOG_INFO("[2/4] linearized: " + std::to_string(unique_chains.size()) + " chains");
    PIRU_PROFILE_STOP(profile, "linearize");

    // -------------------------------------------------------------------------
    // Stage 3: Squigglization + Quantization
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "squigglize");

    piru::signal::FuzzyQuantizerConfig fuzzy_cfg;
    fuzzy_cfg.backend = "rh2";
    auto fuzzy_quantizer = piru::signal::make_fuzzy_quantizer(fuzzy_cfg);
    if (!fuzzy_quantizer) {
        LOG_ERROR("index: failed to create fuzzy quantizer");
        return 1;
    }

    piru::signal::AlignmentQuantizerConfig align_cfg;
    align_cfg.backend = "int16";
    auto alignment_quantizer = piru::signal::make_alignment_quantizer(align_cfg);
    if (!alignment_quantizer) {
        LOG_ERROR("index: failed to create alignment quantizer");
        return 1;
    }

    const auto squiggle_result = piru::index::squigglizeAndQuantize(
        aln_graph, *model, *fuzzy_quantizer, *alignment_quantizer);

    std::size_t total_samples = 0;
    for (const auto& sig : squiggle_result.fuzzy_signals) {
        total_samples += sig.tokens.size();
    }

    LOG_INFO("[3/4] squigglized: " + std::to_string(total_samples) + " signal samples");
    PIRU_PROFILE_STOP(profile, "squigglize");

    // -------------------------------------------------------------------------
    // Stage 4: Seed Extraction & Indexing
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "seeds");

    piru::signal::SeedExtractorConfig extractor_cfg;
    extractor_cfg.backend = "kmer";
    extractor_cfg.k = seed_k;
    extractor_cfg.stride = seed_stride;
    extractor_cfg.qbits = 16;
    auto extractor = piru::signal::make_seed_extractor(extractor_cfg);
    if (!extractor) {
        LOG_ERROR("index: failed to create seed extractor");
        return 1;
    }

    piru::index::SeedBuildConfig seed_cfg;
    seed_cfg.keep_least_frequent_fraction = seed_filter;

    const auto seed_store =
        piru::index::buildSeedStore(squiggle_result.fuzzy_signals, *extractor, seed_cfg);

    LOG_INFO("[4/4] indexed: " + std::to_string(seed_store.size()) + " unique seeds (max_freq=" +
             std::to_string(seed_store.max_hash_frequency()) + ")");
    PIRU_PROFILE_STOP(profile, "seeds");

    // -------------------------------------------------------------------------
    // Stage 5: Serialize Index
    // -------------------------------------------------------------------------

    PIRU_PROFILE_START(profile, "serialize");

    // 1. Create the store objects from the pipeline results.
    piru::index::AdjListGraphStore graph_store(std::move(aln_graph));
    piru::index::VectorSignalStore signal_store(std::move(squiggle_result.alignment_signals));
    // The seed_store is already created.

    // 2. Populate the global metadata.
    piru::io::index::IndexMetadata metadata;
    metadata.graph_flavor = imported.flavor == piru::io::ImportedGraphFlavor::kDbg ? 1 : 2;
    metadata.graph_k = graph_k;
    metadata.pore_k = pore_k;
    metadata.model_name = model->name();
    metadata.fuzzy_quantizer = fuzzy_cfg.backend;
    metadata.align_quantizer = align_cfg.backend;
    metadata.source_path = graph_path;

    // 3. Create the output directory.
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
    }

    // 4. Write the components to disk.
    std::filesystem::path dir_path(output_dir);
    std::string basename = dir_path.filename().string();
    
    piru::io::index::write_graph(dir_path / (basename + ".graph"), graph_store, metadata);
    piru::io::index::write_signals(dir_path / (basename + ".signals"), signal_store, 1.0f, 0.0f);
    piru::io::index::write_seeds(dir_path / (basename + ".seeds"), seed_store);

    LOG_INFO("index written to " + output_dir);
    PIRU_PROFILE_STOP(profile, "serialize");

    PIRU_PROFILE_STOP(profile, "index");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
