// SPDX-License-Identifier: MIT

#include "signal/normalizers/zscore_normalizer.hpp"

#include <cmath>
#include <numeric>
#include <vector>

namespace piru::signal {

ZScoreNormalizer::ZScoreNormalizer(SignalNormalizerConfig config) : config_(std::move(config)) {}

NormalizedSignal ZScoreNormalizer::normalize(const io::RawRead& read, const EventSeries* events) const {
    (void)events;  // Unused for now
    NormalizedSignal normalized;
    auto values = detail::to_picoamps(read);
    if (values.empty()) {
        normalized.sampling_rate_hz = read.sampling_rate_hz;
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
        normalized.samples.push_back(static_cast<float>((static_cast<double>(v) - mean) / stddev));
    }
    normalized.sampling_rate_hz = read.sampling_rate_hz;
    return normalized;
}

const SignalNormalizerConfig& ZScoreNormalizer::config() const { return config_; }

std::string ZScoreNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
