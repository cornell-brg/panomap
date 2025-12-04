#include "commands/index.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "io/graphs/graph.hpp"
#include "io/graphs/graph_loader_factory.hpp"
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
        {'g', "graph", true, "Graph type: dbg (default) or vg"},
        {'m', "model", true, "Pore model (builtin name: r9.4/r10.4 or model file path)"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
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

    // Quick sanity lookup on a canonical homopolymer k-mer for visibility.
    const std::string probe(model->k(), 'A');
    double probe_mean = 0.0;
    const bool probe_ok = model->lookup(probe, probe_mean);

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

    LOG_INFO("graph loaded (" + loader->get_format_name() + "): nodes=" +
             std::to_string(imported.nodes.size()) + ", edges=" +
             std::to_string(imported.edges.size()) + ", paths=" +
             std::to_string(imported.paths.size()));

    if (!imported.nodes.empty()) {
        const auto& n = imported.nodes.front();
        LOG_INFO("first node: id=" + n.id + " len=" + std::to_string(n.sequence.size()));
    }
    if (!imported.edges.empty()) {
        const auto& e = imported.edges.front();
        const auto from_orient = e.from_reverse ? "-" : "+";
        const auto to_orient = e.to_reverse ? "-" : "+";
        std::string overlap_info = e.overlap;
        if (e.overlap_bases.has_value()) {
            overlap_info += " (" + std::to_string(*e.overlap_bases) + "bp)";
        }
        LOG_INFO("first edge: " + e.from + from_orient + " -> " + e.to + to_orient +
                 " overlap=" + overlap_info);
    }
    if (!imported.paths.empty()) {
        const auto& p = imported.paths.front();
        std::string preview;
        const std::size_t limit = std::min<std::size_t>(p.steps.size(), 5);
        for (std::size_t i = 0; i < limit; ++i) {
            const auto& step = p.steps[i];
            preview += (i == 0 ? "" : ",");
            preview += step.segment_id + (step.is_reverse ? "-" : "+");
        }
        if (p.steps.size() > limit) preview += ",...";
        LOG_INFO("first path: name=" + p.name + " steps=" + std::to_string(p.steps.size()) +
                 " [" + preview + "]");
    }

    std::string flavor_str = graph_type;
    if (imported.flavor == piru::io::ImportedGraphFlavor::kUnknown) {
        flavor_str = "unknown(" + graph_type + ")";
    }

    LOG_INFO("index summary (graph=" + flavor_str + ", file=" + graph_path +
             ", model=" + model->name() + ", k=" + std::to_string(model->k()) + ").");
    if (probe_ok) {
        LOG_INFO("model lookup " + probe + " -> mean=" + std::to_string(probe_mean));
    } else {
        LOG_WARN("model lookup failed for probe k-mer: " + probe);
    }

    PIRU_PROFILE_STOP(profile, "index");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
