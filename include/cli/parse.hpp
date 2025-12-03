#pragma once

#include <getopt.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace piru::cli {

struct Option {
    char short_opt;
    std::string long_opt;
    bool requires_arg;
    std::string help;
};

struct ParseConfig {
    std::string usage;
    std::vector<std::string> positional_help;
    std::vector<Option> options;
    std::function<void(const std::string&)> on_error;
};

struct Parsed {
    std::map<std::string, std::string> values;
    std::vector<std::string> positionals;
};

// Build optstring/long options and parse argv using getopt_long.
// Returns false if parsing failed (error already printed via on_error).
bool parse_args(const std::vector<std::string>& args, const ParseConfig& config, Parsed& out);

void print_help(const ParseConfig& config, std::ostream& os);

}  // namespace piru::cli
