#pragma once

#include <cstdlib>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace panomap::cli {

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

// Preset name -> the arguments it stands for.
using PresetTable = std::map<std::string, std::vector<std::string>>;

// Replace `-x <name>` / `-x<name>` / `--preset <name>` / `--preset=<name>` in args
// with the preset's arguments, placed before all other arguments, so that an
// explicit flag later on the command line overrides the preset value.
// Returns false and sets *error for an unknown preset or a missing name.
// If name is given, it receives the last preset named (empty if none).
bool expand_preset(std::vector<std::string>& args, const PresetTable& presets,
                   std::string* error, std::string* name = nullptr);

}  // namespace panomap::cli
