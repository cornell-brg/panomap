#include "commands/index_vg.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

int handle_index_vg(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru index-vg [options]";
    config.options = {
        {'h', "help", false, "Show help"},
        {'r', "reads", true, "Input reads (fastq/bam)"},
        {'g', "graph", true, "Graph description (vg)"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
    };
    config.on_error = [](const std::string&) {
        std::cerr << "index-vg: invalid option\n";
    };

    if (!piru::cli::parse_args(args, config, parsed)) {
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.values.count("help")) {
        piru::cli::print_help(config, std::cout);
        return 0;
    }
    const bool profile = parsed.values.count("profile") > 0;
    PIRU_PROFILE_START(profile, "index-vg");

    LOG_ERROR("index-vg under construction.");
    PIRU_PROFILE_STOP(profile, "index-vg");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
