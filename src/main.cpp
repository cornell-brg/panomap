#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "commands/annotate.hpp"
#include "commands/eval.hpp"
#include "commands/index.hpp"
#include "commands/map.hpp"
#include "commands/mt_test.hpp"
#include "util/metrics.hpp"
#include "util/signal_handlers.hpp"
#include "version.hpp"

struct Command {
    std::string name;
    std::string description;
    int (*handler)(const std::vector<std::string>& args);
};

void print_top_level_usage(const std::vector<Command>& commands) {
    const std::string version = piru_version();
    std::cout << version << "\n\n";
    std::cout << "Usage: piru [--help] [--version] <command> [options]\n";
    std::cout << "Subcommands:\n";
    size_t max_name = 0;
    for (const auto& cmd : commands) {
        if (cmd.name.size() > max_name) {
            max_name = cmd.name.size();
        }
    }
    for (const auto& cmd : commands) {
        std::cout << "  " << cmd.name;
        if (cmd.name.size() < max_name) {
            std::cout << std::string(max_name - cmd.name.size(), ' ');
        }
        std::cout << "  " << cmd.description << "\n";
    }
    std::cout << "\nRun 'piru <command> --help' for command-specific options.\n";
}

int main(int argc, char** argv) {
    piru::install_signal_handlers();

    const double realtime_0 = piru::realtime();
    const double cputime_0 = piru::cputime();

    const std::vector<Command> commands = {
        {"index", "Build index from a pangenome graph.", handle_index},
        {"map", "Map reads against an index (stub).", handle_map},
        {"eval", "Evaluate mapping accuracy against ground truth.", handle_eval},
        {"annotate", "Project BED target regions onto graph node sets.", handle_annotate},
        {"mt-test", "Spawn parallel sleep tasks to test concurrency.", handle_mt_test},
    };

    std::vector<std::string> args(argv + 1, argv + argc);

    int exit_code = 0;

    if (args.empty()) {
        print_top_level_usage(commands);
    } else if (args[0] == "--help" || args[0] == "-h") {
        print_top_level_usage(commands);
    } else if (args[0] == "--version") {
        std::cout << piru_version() << "\n";
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
    std::cerr << "[main] Version: " << PIRU_VERSION << "\n";
    std::cerr << "[main] CMD: ";
    std::cerr << argv[0];
    for (const auto& arg : args) {
        std::cerr << " " << arg;
    }
    std::cerr << "\n[main] Real time: " << std::fixed << std::setprecision(3)
              << (piru::realtime() - realtime_0) << " sec; CPU time: "
              << (piru::cputime() - cputime_0) << " sec; Peak RAM: "
              << (piru::peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB\n\n";

    return exit_code;
}
