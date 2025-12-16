// SPDX-License-Identifier: MIT

#include "signal/normalizers/zscore_normalizer.hpp"

#include <cmath>
#include <numeric>
#include <vector>

namespace piru::signal {

ZScoreNormalizer::ZScoreNormalizer(SignalNormalizerConfig config) : config_(std::move(config)) {}

NormalizedSignal ZScoreNormalizer::normalize(const EventSeries& events) const {
    NormalizedSignal normalized;
    normalized.sampling_rate_hz = events.sampling_rate_hz;

    if (events.events.empty()) {
        return normalized;
    }

    // Collect event means (events already in picoamps)
    std::vector<float> values;
    values.reserve(events.events.size());
    for (const auto& event : events.events) {
        values.push_back(event.mean);
    }

    if (values.empty()) {
        return normalized;
    }

    // Compute global statistics from event means
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
