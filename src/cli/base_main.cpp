/**
 * base_main.cpp
 *
 * Entry point for panomap-base. Wires the base-mode CLI commands and shares
 * the timing/version trailer with panomap.
 *
 * SPDX-License-Identifier: MIT
 */

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "base/commands/index.hpp"
#include "base/commands/inspect.hpp"
#include "base/commands/map.hpp"
#include "cli/parse.hpp"
#include "core/util/metrics.hpp"
#include "core/util/signal_handlers.hpp"
#include "version.hpp"

namespace {

struct Command {
  std::string name;
  std::string description;
  int (*handler)(const std::vector<std::string>& args);
};

void print_top_level_usage(const std::vector<Command>& commands) {
  std::cout << panomap_version() << " (base mode)\n\n";
  std::cout << "Usage: panomap-base [--help] [--version] <command> [options]\n";
  std::cout << "Subcommands:\n";
  std::size_t max_name = 0;
  for (const auto& cmd : commands) {
    if (cmd.name.size() > max_name) max_name = cmd.name.size();
  }
  for (const auto& cmd : commands) {
    std::cout << "  " << cmd.name;
    if (cmd.name.size() < max_name) {
      std::cout << std::string(max_name - cmd.name.size(), ' ');
    }
    std::cout << "  " << cmd.description << "\n";
  }
  std::cout << "\nRun 'panomap-base <command> --help' for command-specific options.\n";
}

}  // namespace

int main(int argc, char** argv) {
  panomap::install_signal_handlers();

  const double realtime_0 = panomap::realtime();
  const double cputime_0 = panomap::cputime();

  const std::vector<Command> commands = {
      {"index", "Build a base-mode minimizer index from a pangenome graph.",
       handle_base_index},
      {"map", "Map basecalled reads against a base-mode index.", handle_base_map},
      {"inspect", "Show base-mode index metadata and stats.", handle_base_inspect},
  };

  std::vector<std::string> args(argv + 1, argv + argc);
  int exit_code = 0;

  if (args.empty() || args[0] == "--help" || args[0] == "-h") {
    print_top_level_usage(commands);
  } else if (args[0] == "--version") {
    std::cout << panomap_version() << "\n";
  } else {
    const std::string& command_name = args[0];
    std::vector<std::string> command_args(args.begin() + 1, args.end());

    bool found = false;
    for (const auto& cmd : commands) {
      if (cmd.name == command_name) {
        exit_code = cmd.handler(command_args);
        found = true;
        break;
      }
    }
    if (!found) {
      std::cerr << "Unknown command: " << command_name << "\n\n";
      print_top_level_usage(commands);
      exit_code = 1;
    }
  }

  std::cerr << "\n--------------------------------\n";
  std::cerr << "[main] Version: " << PANOMAP_VERSION << "\n";
  std::cerr << "[main] CMD: " << argv[0];
  for (const auto& arg : args) {
    std::cerr << " " << arg;
  }
  std::cerr << "\n[main] Real time: " << std::fixed << std::setprecision(3)
            << (panomap::realtime() - realtime_0)
            << " sec; CPU time: " << (panomap::cputime() - cputime_0)
            << " sec; Peak RAM: " << (panomap::peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB\n\n";

  return exit_code;
}
