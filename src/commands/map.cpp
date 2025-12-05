#include "commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "io/reads/read_provider_factory.hpp"
#include "mapping/batch_mapper.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

int handle_map(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru map [options] <reads-path>";
    config.positional_help = {"<reads-path>       Input slow5/blow5 file or directory containing reads"};
    config.options = {
        {'h', "help", false, "Show help"},
        {'t', "threads", true, "Worker threads (-1 = auto)"},
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
    if (parsed.positionals.size() != 1) {
        std::cerr << "map: expected a reads file or directory\n";
        piru::cli::print_help(config, std::cerr);
        return 1;
    }

    const bool profile = parsed.values.count("profile") > 0;
    const int num_threads = [&]() {
        auto it = parsed.values.find("threads");
        if (it == parsed.values.end()) return -1;
        try {
            return std::stoi(it->second);
        } catch (...) {
            LOG_WARN("map: invalid --threads value '" + it->second + "', using auto");
            return -1;
        }
    }();
    PIRU_PROFILE_START(profile, "map");

    const std::string reads_path = parsed.positionals[0];
    std::vector<std::filesystem::path> files;
    std::error_code ec;
    if (std::filesystem::is_directory(reads_path, ec)) {
        for (const auto& entry : std::filesystem::recursive_directory_iterator(reads_path)) {
            if (!entry.is_regular_file()) continue;
            const auto ext = entry.path().extension().string();
            const std::string ext_lower = [&]() {
                std::string s = ext;
                std::transform(s.begin(), s.end(), s.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                return s;
            }();
            if (ext_lower == ".blow5" || ext_lower == ".slow5") {
                files.push_back(entry.path());
            }
        }
    } else {
        files.push_back(reads_path);
    }

    if (files.empty()) {
        LOG_ERROR("map: no supported read files found under '" + reads_path +
                  "' (expects .blow5/.slow5)");
        return 1;
    }

    piru::mapping::BatchMapperConfig map_config;
    map_config.num_threads = num_threads;

    std::size_t total_reads = 0;
    std::size_t total_batches = 0;
    std::size_t files_processed = 0;

    for (const auto& f : files) {
        auto provider = piru::io::make_read_provider(f.string());
        if (!provider) {
            LOG_WARN("map: unsupported read format for '" + f.string() + "', skipping");
            continue;
        }
        ++files_processed;
        LOG_INFO("map: processing '" + f.string() + "' (format=" + provider->get_format_name() +
                 ").");

        piru::mapping::BatchMapper mapper(*provider, map_config, std::cout);
        const auto stats = mapper.process_all();
        total_reads += stats.reads_processed;
        total_batches += stats.batches;
    }

    if (files_processed == 0) {
        LOG_ERROR("map: no readable input files");
        return 1;
    }

    LOG_INFO("map: done. files=" + std::to_string(files_processed) +
             ", batches=" + std::to_string(total_batches) +
             ", reads=" + std::to_string(total_reads));

    PIRU_PROFILE_STOP(profile, "map");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
