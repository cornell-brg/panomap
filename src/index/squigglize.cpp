// SPDX-License-Identifier: MIT
#include "index/squigglize.hpp"

#include <cmath>
#include <limits>

#include "util/logging.hpp"

namespace piru::index {
namespace {

// Check if k-mer contains N or other non-ACGT bases
bool hasNBase(const std::string_view& kmer) {
    for (char c : kmer) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return true;
        }
    }
    return false;
}

float safeNormalize(float value, double mean, double stddev) {
    // If all values are identical, normalize to 0
    if (stddev < 1e-6) return 0.0f;

    // Z-score normalization with outlier clipping [-3, 3]
    float z_score = static_cast<float>((value - mean) / stddev);
    if (z_score < -3.0f) z_score = -3.0f;
    if (z_score > 3.0f) z_score = 3.0f;
    return z_score;
}

}  // namespace

SquiggleResult squigglizeAndQuantize(const AlnGraph& graph,
                                     const piru::io::KmerModel& model,
                                     const piru::signal::FuzzyQuantizer& fuzzy_quantizer,
                                     const piru::signal::AlignmentQuantizer& alignment_quantizer) {
    const int k = model.k();
    SquiggleResult result;
    result.fuzzy_signals.resize(graph.nodeCount());
    result.alignment_signals.resize(graph.nodeCount());

    std::vector<std::vector<float>> raw_signals(graph.nodeCount());
    double sum = 0.0;
    std::size_t count = 0;

    // Generate expected signals per node and accumulate global stats (first pass: compute mean).
    for (std::size_t node_id = 0; node_id < graph.nodeCount(); ++node_id) {
        const auto& node = graph.node(node_id);
        const std::string& seq = node.sequence;
        if (seq.size() < static_cast<std::size_t>(k)) {
            continue;
        }

        auto& samples = raw_signals[node_id];
        samples.reserve(seq.size() - k + 1);
        std::size_t n_kmer_count = 0;
        for (std::size_t i = 0; i + static_cast<std::size_t>(k) <= seq.size(); ++i) {
            const std::string_view kmer(seq.data() + i, static_cast<std::size_t>(k));

            // Check for N bases or other non-ACGT characters
            if (hasNBase(kmer)) {
                // Emit NaN sentinel for N-containing k-mers
                samples.push_back(std::numeric_limits<float>::quiet_NaN());
                ++n_kmer_count;
                continue;
            }

            double mean = 0.0;
            if (!model.lookup(std::string(kmer), mean)) {
                LOG_WARN("Missing k-mer in model at node " + node.label + " pos " +
                         std::to_string(i) + ": " + std::string(kmer));
                mean = 0.0;
            }
            const float val = static_cast<float>(mean);
            samples.push_back(val);
            sum += val;
            ++count;
        }

        // // Log summary if node has N bases
        // if (n_kmer_count > 0) {
        //     LOG_INFO("Node " + node.label + ": " + std::to_string(n_kmer_count) +
        //              " k-mers with N bases (marked as NaN)");
        // }
    }

    const double global_mean = (count == 0) ? 0.0 : sum / static_cast<double>(count);

    // Second pass: compute variance using two-pass algorithm for numerical stability.
    // Skip NaN values when computing variance.
    double variance = 0.0;
    if (count > 0) {
        for (const auto& samples : raw_signals) {
            for (const auto val : samples) {
                if (std::isnan(val)) continue;  // Skip NaN sentinels
                const double diff = static_cast<double>(val) - global_mean;
                variance += diff * diff;
            }
        }
        variance /= static_cast<double>(count);
    }
    const double global_std = (variance > 0.0) ? std::sqrt(variance) : 0.0;

    // Normalize and quantize.
    for (std::size_t node_id = 0; node_id < graph.nodeCount(); ++node_id) {
        const auto& samples = raw_signals[node_id];
        piru::signal::NormalizedSignal normalized;
        normalized.sampling_rate_hz = 0.0f;
        normalized.samples.reserve(samples.size());
        for (const auto val : samples) {
            // Pass through NaN unchanged; normalize valid values
            if (std::isnan(val)) {
                normalized.samples.push_back(val);
            } else {
                normalized.samples.push_back(safeNormalize(val, global_mean, global_std));
            }
        }

        result.fuzzy_signals[node_id] = fuzzy_quantizer.quantize(normalized, nullptr);
        result.alignment_signals[node_id] = alignment_quantizer.quantize(normalized, nullptr);
    }

    result.raw_signals = std::move(raw_signals);
    return result;
}

}  // namespace piru::index
