// SPDX-License-Identifier: MIT

#include "index/path_walk_indexer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "util/logging.hpp"

namespace piru::index {

namespace {

/* Check if k-mer contains N or other non-ACGT bases. */
bool hasNBase(const std::string_view& kmer) {
    for (char c : kmer) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return true;
        }
    }
    return false;
}

/* Z-score normalization with outlier clipping to [-3, 3]. */
float safeNormalize(float value, double mean, double stddev) {
    if (stddev < 1e-6) return 0.0f;
    float z_score = static_cast<float>((value - mean) / stddev);
    return std::clamp(z_score, -3.0f, 3.0f);
}

/* Track where each node starts in the concatenated path sequence. */
struct NodeBoundary {
    std::size_t node_id;
    std::size_t base_start;
};

/* Map signal position back to (node_id, local_offset). */
std::pair<std::size_t, std::size_t> signalPositionToNodeOffset(
    std::size_t signal_pos,
    const std::vector<NodeBoundary>& boundaries) {
    for (auto it = boundaries.rbegin(); it != boundaries.rend(); ++it) {
        if (signal_pos >= it->base_start) {
            return {it->node_id, signal_pos - it->base_start};
        }
    }
    return {boundaries[0].node_id, signal_pos};
}

}  // namespace

PathWalkIndexResult pathWalkIndex(
    const AlnGraph& graph,
    const io::KmerModel& model,
    const signal::FuzzyQuantizer& fuzzy_quantizer,
    const signal::SeedExtractor& extractor,
    const PathWalkIndexConfig& config)
{
    PathWalkIndexResult result;
    result.seed_store = std::make_unique<HashSeedStore>();
    result.linearization_coords.resize(graph.nodeCount());
    result.path_lengths.resize(graph.pathCount(), 0);

    const int pore_k = model.k();

    // Store extractor metadata for serialization
    result.seed_store->set_extractor_name(extractor.name());
    std::map<std::string, std::string> params;
    const auto& cfg = extractor.config();
    params["backend"] = cfg.backend;
    params["k"] = std::to_string(cfg.k);
    params["stride"] = std::to_string(cfg.stride);
    params["qbits"] = std::to_string(cfg.qbits);
    result.seed_store->set_params(std::move(params));
    result.seed_store->set_filter_fraction(config.seed_filter);

    // Process each path independently
    for (std::size_t path_idx = 0; path_idx < graph.pathCount(); ++path_idx) {
        const auto& path = graph.paths()[path_idx];

        // Step 1: Concatenate node sequences, track boundaries
        std::string path_sequence;
        std::vector<NodeBoundary> boundaries;

        for (const auto& step : path.steps) {
            std::size_t node_id = std::stoull(step.node_id);
            if (node_id >= graph.nodeCount()) continue;

            const auto& node = graph.node(node_id);
            boundaries.push_back({node_id, path_sequence.size()});
            path_sequence += node.sequence;
        }

        if (path_sequence.size() < static_cast<std::size_t>(pore_k)) {
            LOG_WARN("Path " + path.name + " too short for k=" + std::to_string(pore_k) + ", skipping...");
            continue;
        }

        // Step 2: Squigglize (slide k-mer window, k-mers can cross node boundaries)
        std::vector<float> raw_signal;
        raw_signal.reserve(path_sequence.size() - pore_k + 1);

        double sum = 0.0;
        std::size_t count = 0;
        std::string kmer_buf(static_cast<std::size_t>(pore_k), '\0');

        for (std::size_t i = 0; i + static_cast<std::size_t>(pore_k) <= path_sequence.size(); ++i) {
            std::copy_n(path_sequence.data() + i, static_cast<std::size_t>(pore_k), kmer_buf.begin());

            if (hasNBase(kmer_buf)) {
                raw_signal.push_back(std::numeric_limits<float>::quiet_NaN());
                continue;
            }

            double mean = 0.0;
            if (!model.lookup(kmer_buf, mean)) {
                mean = 0.0;
            }
            float val = static_cast<float>(mean);
            raw_signal.push_back(val);
            sum += val;
            ++count;
        }

        if (count == 0) continue;

        // Step 3: Compute per-path normalization stats
        const double path_mean = sum / static_cast<double>(count);

        double variance = 0.0;
        for (float val : raw_signal) {
            if (!std::isnan(val)) {
                double diff = static_cast<double>(val) - path_mean;
                variance += diff * diff;
            }
        }
        variance /= static_cast<double>(count);
        const double path_std = (variance > 0.0) ? std::sqrt(variance) : 0.0;

        // Step 4: Normalize and fuzzy quantize
        signal::NormalizedSignal normalized;
        normalized.samples.reserve(raw_signal.size());
        for (float val : raw_signal) {
            if (std::isnan(val)) {
                normalized.samples.push_back(val);
            } else {
                normalized.samples.push_back(safeNormalize(val, path_mean, path_std));
            }
        }

        auto fuzzy = fuzzy_quantizer.quantize(normalized);
        result.path_lengths[path_idx] = fuzzy.tokens.size();
        result.total_path_length += fuzzy.tokens.size();

        // Step 5: Record linear coordinates (node_id -> path positions for chaining)
        for (const auto& boundary : boundaries) {
            result.linearization_coords[boundary.node_id].emplace_back(
                path_idx, static_cast<std::int64_t>(boundary.base_start));
        }

        // Step 6: Extract seeds and insert into index
        const auto seeds = extractor.extract(fuzzy);
        result.seeds_extracted += seeds.seeds.size();

        for (const auto& seed : seeds.seeds) {
            auto [node_id, local_offset] = signalPositionToNodeOffset(seed.position, boundaries);
            result.seed_store->insert(seed.hash, SeedHit{node_id, local_offset, seed.length});
        }
    }

    // Dedup seeds from shared regions across paths
    result.seed_store->deduplicate();
    result.seeds_unique = result.seed_store->size();

    // Compute max frequency
    std::size_t max_freq = 0;
    for (const auto& [hash, hits] : result.seed_store->data()) {
        max_freq = std::max(max_freq, hits.size());
    }
    result.seed_store->set_max_hash_frequency(max_freq);

    // Compute frequency threshold (seeds NOT deleted, threshold used at query time)
    std::vector<std::size_t> frequencies;
    frequencies.reserve(result.seed_store->size());
    for (const auto& [hash, hits] : result.seed_store->data()) {
        frequencies.push_back(hits.size());
    }

    std::size_t threshold = max_freq + 1;
    if (!frequencies.empty() && config.seed_filter < 1.0) {
        std::sort(frequencies.begin(), frequencies.end());
        double fraction = std::clamp(config.seed_filter, 0.0, 1.0);
        std::size_t pos = static_cast<std::size_t>(frequencies.size() * fraction);
        pos = std::min(pos, frequencies.size() - 1);
        threshold = frequencies[pos] + 1;
    }
    result.seed_store->set_frequency_threshold(threshold);

    LOG_INFO("PathWalkIndex: " + std::to_string(result.seeds_extracted) +
             " seeds extracted, " + std::to_string(result.seeds_unique) +
             " unique (max_freq=" + std::to_string(max_freq) + ")");

    return result;
}

}  // namespace piru::index
