#include "commands/map.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "io/reads/read_provider_factory.hpp"
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
        {'t', "threads", true, "Worker threads (unused in demo)"},
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

    // Choose the first file and its provider as the runtime format for this invocation.
    auto provider = piru::io::make_read_provider(files.front().string());
    if (!provider) {
        LOG_ERROR("map: unsupported read format for '" + files.front().string() + "'");
        return 1;
    }
    const std::string chosen_ext = files.front().extension().string();
    std::vector<std::filesystem::path> filtered;
    for (const auto& f : files) {
        if (f.extension().string() == chosen_ext) {
            filtered.push_back(f);
        }
    }

    std::size_t total_reads = 0;
    for (const auto& f : filtered) {
        // Recreate provider per file to keep the API simple for now.
        auto p = piru::io::make_read_provider(f.string());
        if (!p) {
            LOG_ERROR("map: unsupported read format for '" + f.string() + "'");
            continue;
        }
        std::size_t count = 0;
        piru::io::RawRead read;
        while (p->get_next(read)) {
            ++count;
            std::cout << read.read_id << " len=" << read.len_raw_signal << " range=" << read.range
                      << " digitisation=" << read.digitisation << " offset=" << read.offset
                      << " samplerate=" << read.sampling_rate_hz << "Hz\n";
        }
        total_reads += count;
        LOG_INFO("map demo listing '" + f.string() + "' (format=" + p->get_format_name() +
                 ", reads=" + std::to_string(count) + ").");
    }

    LOG_INFO("map demo total reads listed: " + std::to_string(total_reads) +
             " (using format " + provider->get_format_name() + ").");

    PIRU_PROFILE_STOP(profile, "map");
    if (profile) piru::timing::report(std::cerr);
    return 0;
}
