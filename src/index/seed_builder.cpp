// SPDX-License-Identifier: MIT
#include "index/seed_builder.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace piru::index {

HashSeedStore buildSeedStore(const std::vector<piru::signal::FuzzyQuantizedSignal>& signals,
                             const piru::signal::SeedExtractor& extractor,
                             const SeedBuildConfig& config) {
    HashSeedStore store;

    // Record extractor metadata for compatibility checks during mapping.
    store.set_extractor_name(extractor.name());
    std::map<std::string, std::string> params;
    const auto& cfg = extractor.config();
    params["backend"] = cfg.backend;
    params["k"] = std::to_string(cfg.k);
    params["stride"] = std::to_string(cfg.stride);
    params["qbits"] = std::to_string(cfg.qbits);
    if (cfg.window > 0) {
        params["window"] = std::to_string(cfg.window);
    }
    if (!cfg.params.empty()) {
        params["params"] = cfg.params;
    }
    store.set_params(std::move(params));
    store.set_filter_fraction(config.keep_least_frequent_fraction);

    // Step 1: Extract seeds from all nodes and populate hash table
    for (std::size_t node_id = 0; node_id < signals.size(); ++node_id) {
        const auto& signal = signals[node_id];
        const auto seeds = extractor.extract(signal, nullptr);
        for (const auto& seed : seeds.seeds) {
            store.insert(seed.hash, SeedHit{node_id, seed.position});
        }
    }

    // Step 2: Compute max hash frequency
    std::size_t max_freq = 0;
    for (const auto& [hash, hits] : store.data()) {
        if (hits.size() > max_freq) {
            max_freq = hits.size();
        }
    }
    store.set_max_hash_frequency(max_freq);

    // Step 3: Compute frequency threshold (GraphAligner-style percentile filtering)
    std::vector<std::size_t> frequencies;
    frequencies.reserve(store.size());
    for (const auto& [hash, hits] : store.data()) {
        frequencies.push_back(hits.size());
    }

    std::size_t threshold = 0;
    if (!frequencies.empty() && config.keep_least_frequent_fraction < 1.0) {
        std::sort(frequencies.begin(), frequencies.end());
        double fraction = config.keep_least_frequent_fraction;
        if (fraction < 0.0) fraction = 0.0;
        if (fraction > 1.0) fraction = 1.0;

        std::size_t pos = static_cast<std::size_t>(frequencies.size() * fraction);
        if (pos >= frequencies.size()) {
            pos = frequencies.size() - 1;
        }
        // Threshold is percentile + 1 (GraphAligner convention: keep hashes with freq <= threshold)
        threshold = frequencies[pos] + 1;
    } else {
        // Keep all seeds: threshold = max_freq + 1
        threshold = max_freq + 1;
    }
    store.set_frequency_threshold(threshold);

    return store;
}

}  // namespace piru::index
