// SPDX-License-Identifier: MIT
// Interface for normalizing raw read signals.

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "io/reads/read_provider.hpp"
#include "signal/signal_types.hpp"

namespace piru::signal {

struct SignalNormalizerConfig {
    std::string backend{"zscore"};
};

class SignalNormalizer {
public:
    virtual ~SignalNormalizer() = default;

    virtual NormalizedSignal normalize(const io::RawRead& read, const EventSeries* events) const = 0;
    virtual const SignalNormalizerConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using SignalNormalizerPtr = std::unique_ptr<SignalNormalizer>;

namespace detail {

// Convert raw ADC samples to picoamps using read metadata.
inline std::vector<float> to_picoamps(const io::RawRead& read) {
    const float raw_unit = (read.digitisation == 0.0f) ? 1.0f : (read.range / read.digitisation);
    std::vector<float> pa;
    pa.reserve(read.raw_signal.size());
    for (auto value : read.raw_signal) {
        pa.push_back((static_cast<float>(value) + read.offset) * raw_unit);
    }
    return pa;
}

// In-place median helper (modifies input buffer).
inline float median_inplace(std::vector<float>& vec) {
    const std::size_t n = vec.size();
    if (n == 0) return 0.0f;
    const std::size_t half = n / 2;
    std::nth_element(vec.begin(), vec.begin() + half, vec.end());
    float med = vec[half];
    if ((n % 2) == 0) {
        std::nth_element(vec.begin(), vec.begin() + half - 1, vec.end());
        med = 0.5f * (med + vec[half - 1]);
    }
    return med;
}

}  // namespace detail

}  // namespace piru::signal
