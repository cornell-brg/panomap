// SPDX-License-Identifier: MIT
#include "index/seed_builder.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "util/logging.hpp"

namespace piru::index {

namespace {

// Build label → numeric node id map for path-guided seeding.
std::unordered_map<std::string, std::size_t> buildLabelToIdMap(const AlnGraph& graph) {
    std::unordered_map<std::string, std::size_t> label_to_id;
    for (std::size_t i = 0; i < graph.nodeCount(); ++i) {
        const auto& node = graph.node(i);
        label_to_id[node.label] = node.id;
    }
    return label_to_id;
}

// Map a position in concatenated path signal back to (node_id, local_offset).
// node_offsets: [(node_id, cumulative_start), ...] sorted by cumulative_start.
std::pair<std::size_t, std::size_t> positionToNodeOffset(
    std::size_t position,
    const std::vector<std::pair<std::size_t, std::size_t>>& node_offsets) {
    // Binary search for the node containing this position.
    // Find the last node whose cumulative_start <= position.
    auto it = std::upper_bound(
        node_offsets.begin(), node_offsets.end(), position,
        [](std::size_t pos, const std::pair<std::size_t, std::size_t>& entry) {
            return pos < entry.second;
        });
    if (it != node_offsets.begin()) {
        --it;
    }
    std::size_t node_id = it->first;
    std::size_t local_offset = position - it->second;
    return {node_id, local_offset};
}

}  // namespace

HashSeedStore buildSeedStore(const AlnGraph* graph,
                             const std::vector<piru::signal::FuzzyQuantizedSignal>& signals,
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

    // Step 1: Extract seeds and populate hash table
    std::size_t total_seeds_extracted = 0;

    if (graph == nullptr) {
        // Node-based seeding: extract from each node independently
        for (std::size_t node_id = 0; node_id < signals.size(); ++node_id) {
            const auto& signal = signals[node_id];
            const auto seeds = extractor.extract(signal);
            total_seeds_extracted += seeds.seeds.size();
            for (const auto& seed : seeds.seeds) {
                store.insert(seed.hash, SeedHit{node_id, seed.position, seed.length});
            }
        }

        LOG_INFO("Seed extraction (node-based): extracted=" +
                 std::to_string(total_seeds_extracted) + " seeds from " +
                 std::to_string(signals.size()) + " nodes, unique=" +
                 std::to_string(store.size()) + " (" +
                 std::to_string(static_cast<int>(
                     total_seeds_extracted > 0
                         ? (100.0 * store.size() / total_seeds_extracted)
                         : 0.0)) +
                 "%)");
    } else {
        // Path-guided seeding: walk paths, extract seeds that can cross node boundaries
        const auto label_to_id = buildLabelToIdMap(*graph);

        for (const auto& path : graph->paths()) {
            // Build node offset table and concatenate signals for this path
            std::vector<std::pair<std::size_t, std::size_t>> node_offsets;
            std::vector<std::int16_t> path_tokens;

            for (const auto& step : path.steps) {
                auto it = label_to_id.find(step.node_id);
                if (it == label_to_id.end()) {
                    continue;  // Node not found, skip
                }
                std::size_t node_id = it->second;
                if (node_id >= signals.size()) {
                    continue;  // Signal not available, skip
                }

                node_offsets.emplace_back(node_id, path_tokens.size());
                const auto& signal = signals[node_id];
                path_tokens.insert(path_tokens.end(), signal.tokens.begin(),
                                   signal.tokens.end());
            }

            if (path_tokens.empty()) {
                continue;
            }

            // Extract seeds from concatenated path signal
            piru::signal::FuzzyQuantizedSignal path_signal{std::move(path_tokens)};
            const auto seeds = extractor.extract(path_signal);
            total_seeds_extracted += seeds.seeds.size();

            // Map seed positions back to (node_id, local_offset) and insert
            for (const auto& seed : seeds.seeds) {
                auto [node_id, local_offset] =
                    positionToNodeOffset(seed.position, node_offsets);
                store.insert(seed.hash, SeedHit{node_id, local_offset, seed.length});
            }
        }

        // Deduplicate seeds from shared regions across paths
        store.deduplicate();

        LOG_INFO("Seed extraction (path-guided): extracted=" +
                 std::to_string(total_seeds_extracted) + " seeds from " +
                 std::to_string(graph->pathCount()) + " paths, unique=" +
                 std::to_string(store.size()) + " (after dedup)");
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
