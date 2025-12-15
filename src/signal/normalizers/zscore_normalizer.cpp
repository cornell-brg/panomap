// SPDX-License-Identifier: MIT

#include "signal/normalizers/zscore_normalizer.hpp"

#include <cmath>
#include <numeric>
#include <vector>

namespace piru::signal {

ZScoreNormalizer::ZScoreNormalizer(SignalNormalizerConfig config) : config_(std::move(config)) {}

NormalizedSignal ZScoreNormalizer::normalize(const io::RawRead& read, const EventSeries* events) const {
    NormalizedSignal normalized;
    normalized.sampling_rate_hz = read.sampling_rate_hz;

    // If events are not provided, fall back to raw signal conversion
    if (!events || events->events.empty()) {
        auto values = detail::to_picoamps(read);
        if (values.empty()) {
            return normalized;
        }

        const double sum = std::accumulate(values.begin(), values.end(), 0.0);
        const double mean = sum / static_cast<double>(values.size());

        double var = 0.0;
        for (const auto v : values) {
            const double diff = static_cast<double>(v) - mean;
            var += diff * diff;
        }
        var /= static_cast<double>(values.size());
        const double stddev = (var > 0.0) ? std::sqrt(var) : 1.0;

        normalized.samples.reserve(values.size());
        for (const auto v : values) {
            float z_score = static_cast<float>((static_cast<double>(v) - mean) / stddev);

            // Apply outlier clipping if enabled
            if (config_.clip_outliers) {
                if (z_score < config_.clip_min) z_score = config_.clip_min;
                if (z_score > config_.clip_max) z_score = config_.clip_max;
            }

            normalized.samples.push_back(z_score);
        }
        return normalized;
    }

    // Reconstruct signal from events (events already in picoamps)
    std::vector<float> values;
    values.reserve(read.len_raw_signal);

    for (const auto& event : events->events) {
        // Fill event duration with event mean value
        for (std::size_t i = 0; i < event.length; ++i) {
            values.push_back(event.mean);
        }
    }

    if (values.empty()) {
        return normalized;
    }

    // Compute global statistics from reconstructed signal
    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    const double mean = sum / static_cast<double>(values.size());

    double var = 0.0;
    for (const auto v : values) {
        const double diff = static_cast<double>(v) - mean;
        var += diff * diff;
    }
    var /= static_cast<double>(values.size());
    const double stddev = (var > 0.0) ? std::sqrt(var) : 1.0;

    // Z-score normalize
    normalized.samples.reserve(values.size());
    for (const auto v : values) {
        float z_score = static_cast<float>((static_cast<double>(v) - mean) / stddev);

        // Apply outlier clipping if enabled
        if (config_.clip_outliers) {
            if (z_score < config_.clip_min) z_score = config_.clip_min;
            if (z_score > config_.clip_max) z_score = config_.clip_max;
        }

        normalized.samples.push_back(z_score);
    }

    return normalized;
}

const SignalNormalizerConfig& ZScoreNormalizer::config() const { return config_; }

std::string ZScoreNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
