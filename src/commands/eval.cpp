#include "commands/eval.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "util/logging.hpp"
#include "version.hpp"

int handle_eval(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru eval [options]";
    config.options = {
        {'h', "help", false, "Show help"},
        {'t', "truth", true, "Ground truth/reference annotations"},
        {'c', "calls", true, "Calls to evaluate for enrichment/depletion"},
    };
    config.on_error = [](const std::string&) { std::cerr << "eval: invalid option\n"; };

    if (!piru::cli::parse_args(args, config, parsed)) {
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.values.count("help")) {
        piru::cli::print_help(config, std::cout);
        return 0;
    }
    LOG_ERROR("eval under construction.");
    return 0;
}
