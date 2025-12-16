// SPDX-License-Identifier: MIT

#include "signal/normalizers/median_mad_normalizer.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace piru::signal {

namespace {

constexpr float kMadScale = 1.4826f;

}  // namespace

MedianMadNormalizer::MedianMadNormalizer(SignalNormalizerConfig config)
    : config_(std::move(config)) {}

NormalizedSignal MedianMadNormalizer::normalize(const EventSeries& events) const {
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

    std::vector<float> scratch = values;
    const float med = detail::median_inplace(scratch);
    for (auto& v : scratch) {
        v = std::abs(v - med);
    }
    float mad = detail::median_inplace(scratch) * kMadScale;
    if (mad == 0.0f) {
        mad = 1.0f;
    }

    normalized.samples.reserve(values.size());
    for (const auto v : values) {
        float norm_value = (v - med) / mad;

        // Apply outlier clipping if enabled
        if (config_.clip_outliers) {
            if (norm_value < config_.clip_min) norm_value = config_.clip_min;
            if (norm_value > config_.clip_max) norm_value = config_.clip_max;
        }

        normalized.samples.push_back(norm_value);
    }
    return normalized;
}

const SignalNormalizerConfig& MedianMadNormalizer::config() const { return config_; }

std::string MedianMadNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
