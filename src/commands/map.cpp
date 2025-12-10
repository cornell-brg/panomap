#include "commands/map.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "io/index/serialization.hpp"
#include "io/reads/read_provider_factory.hpp"
#include "mapping/batch_mapper.hpp"
#include "util/logging.hpp"
#include "util/timing.hpp"
#include "version.hpp"

namespace {

piru::signal::SeedExtractorConfig seed_config_from_store(const piru::index::HashSeedStore& store) {
    piru::signal::SeedExtractorConfig cfg;
    cfg.backend = store.extractor_name();
    const auto& params = store.params();
    auto get_u64 = [&params](const std::string& key, std::size_t default_val) -> std::size_t {
        auto it = params.find(key);
        if (it == params.end()) return default_val;
        return std::stoull(it->second);
    };
    cfg.k = get_u64("k", cfg.k);
    cfg.stride = get_u64("stride", cfg.stride);
    cfg.window = get_u64("window", cfg.window);
    auto it_qbits = params.find("qbits");
    if (it_qbits != params.end()) cfg.qbits = static_cast<std::uint32_t>(std::stoul(it_qbits->second));
    auto it_params = params.find("params");
    if (it_params != params.end()) cfg.params = it_params->second;
    return cfg;
}

void apply_index_compatibility(const piru::io::index::LoadedIndex& loaded,
                               piru::mapping::BatchMapperConfig& map_config) {
    // Fuzzy quantizer backend must match index metadata.
    if (!loaded.metadata.fuzzy_quantizer.empty()) {
        map_config.fuzzy_config.backend = loaded.metadata.fuzzy_quantizer;
    }

    // Seed extractor config must match the SeedStore metadata.
    const auto* seed_store = dynamic_cast<piru::index::HashSeedStore*>(loaded.seeds.get());
    if (!seed_store) {
        throw std::runtime_error("Unsupported SeedStore backend in index (expected hash).");
    }
    auto index_seed_cfg = seed_config_from_store(*seed_store);
    if (index_seed_cfg.backend.empty()) {
        throw std::runtime_error("SeedStore extractor name missing in index.");
    }

    // If map config was left at defaults, adopt index config; otherwise require an exact match.
    if (map_config.seed_config != index_seed_cfg) {
        // Adopt automatically only when map config appears default-ish.
        if (map_config.seed_config.backend == "kmer" && map_config.seed_config.k == 10 &&
            map_config.seed_config.stride == 1 && map_config.seed_config.qbits == 4 &&
            map_config.seed_config.window == 0 && map_config.seed_config.params.empty()) {
            LOG_INFO("Using seed extractor settings from index: backend=" + index_seed_cfg.backend +
                     ", k=" + std::to_string(index_seed_cfg.k) +
                     ", stride=" + std::to_string(index_seed_cfg.stride) +
                     ", qbits=" + std::to_string(index_seed_cfg.qbits));
            map_config.seed_config = index_seed_cfg;
        } else {
            throw std::runtime_error(
                "Seed extractor config mismatch between map settings and index (expected backend '" +
                index_seed_cfg.backend + "').");
        }
    }
}

}  // namespace

int handle_map(const std::vector<std::string>& args) {
    piru::cli::Parsed parsed;
    piru::cli::ParseConfig config;
    config.usage = "Usage: piru map [options] --index <index-dir> <reads-path>";
    config.positional_help = {"<reads-path>       Input slow5/blow5 file or directory containing reads"};
    config.options = {
        {'h', "help", false, "Show help"},
        {'i', "index", true, "Path to index directory"},
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
    if (!parsed.values.count("index")) {
        LOG_ERROR("map: missing required --index <index-dir>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }
    if (parsed.positionals.size() != 1) {
        LOG_ERROR("map: missing required <reads-path>");
        piru::cli::print_help(config, std::cerr);
        return 1;
    }

    const std::string index_path = parsed.values.at("index");
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

    LOG_INFO("loading index from " + index_path);
    auto loaded_index = piru::io::index::load_index(index_path);
    if (!loaded_index.graph) {
        LOG_ERROR("map: failed to load index from '" + index_path + "'");
        return 1;
    }
    if (!loaded_index.seeds) {
        LOG_ERROR("map: index is missing seeds; cannot perform lookups");
        return 1;
    }
    LOG_INFO("index loaded: " + std::to_string(loaded_index.graph->nodeCount()) + " nodes");
    LOG_INFO("index metadata: fuzzy=" + loaded_index.metadata.fuzzy_quantizer +
             ", align=" + loaded_index.metadata.align_quantizer +
             ", seeds=" + loaded_index.seeds->extractor_name());

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
    map_config.seed_store = loaded_index.seeds.get();
    map_config.graph_store = loaded_index.graph.get();

    // Align map-side settings with index metadata before processing.
    apply_index_compatibility(loaded_index, map_config);

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
