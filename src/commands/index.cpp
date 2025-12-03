#include "commands/index.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
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

    LOG_ERROR("index under construction (graph=" + graph_type + ", file=" + graph_path +
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
