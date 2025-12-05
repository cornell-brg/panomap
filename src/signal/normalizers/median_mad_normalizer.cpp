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

NormalizedSignal MedianMadNormalizer::normalize(const io::RawRead& read, const EventSeries* events) const {
    (void)events;  // Unused for now
    NormalizedSignal normalized;
    auto values = detail::to_picoamps(read);
    if (values.empty()) {
        normalized.sampling_rate_hz = read.sampling_rate_hz;
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
        normalized.samples.push_back((v - med) / mad);
    }
    normalized.sampling_rate_hz = read.sampling_rate_hz;
    return normalized;
}

const SignalNormalizerConfig& MedianMadNormalizer::config() const { return config_; }

std::string MedianMadNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
