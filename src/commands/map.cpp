#include "commands/map.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

int handle_map(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru map [options]";
    config.options = {
        {'h', "help", false, "Show help"},
        {'r', "reads", true, "Reads to map"},
        {'i', "index", true, "Index to query"},
        {'t', "threads", true, "Worker threads"},
        {'p', "profile", false, "Emit timing profile (tree)"},
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
    const bool profile = parsed.values.count("profile") > 0;
    PIRU_PROFILE_START(profile, "map");

    LOG_ERROR("map under construction.");
    PIRU_PROFILE_STOP(profile, "map");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
