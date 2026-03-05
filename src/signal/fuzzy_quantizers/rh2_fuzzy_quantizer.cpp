// SPDX-License-Identifier: MIT

#include "signal/fuzzy_quantizers/rh2_fuzzy_quantizer.hpp"

#include <algorithm>
#include <cmath>

namespace piru::signal {

namespace {

std::uint32_t dynamic_quantize(float value, const FuzzyQuantizerConfig& cfg) {
    const float fine_min = cfg.fine_min;
    const float fine_max = cfg.fine_max;
    const float fine_range = cfg.fine_range;
    const std::uint32_t n_buckets = 1u << cfg.qbits;

    if (n_buckets == 0u) return 0u;

    float min_val = -3.0f;
    float max_val = 3.0f;
    float range = max_val - min_val;
    float coarse_coef1 = (1.0f - fine_range) * 0.5f;
    float coarse_coef2 = fine_range + coarse_coef1;

    float signal = std::clamp(value, min_val, max_val);

    // Min-max normalization to [0, 1]
    float normalized = (signal - min_val) / range;

    // Fine range mapping
    float a = (fine_min - min_val) / range;
    float b = (fine_max - min_val) / range;

    // Conditional quantization
    float quantized = fine_max;
    if (signal >= fine_min && signal <= fine_max) {
        quantized = fine_range * ((normalized - a) / (b - a));
    } else {
        if (normalized < 0.5f) {
            quantized = fine_range + coarse_coef1 * normalized;
        } else {
            quantized = coarse_coef2 + coarse_coef1 * normalized;
        }
    }

    std::uint32_t qnt_val = static_cast<std::uint32_t>(quantized * (n_buckets - 1));
    return qnt_val;
}

}  // namespace

Rh2FuzzyQuantizer::Rh2FuzzyQuantizer(FuzzyQuantizerConfig config) : config_(std::move(config)) {}

FuzzyQuantizedSignal Rh2FuzzyQuantizer::quantize(const NormalizedSignal& signal) const {
    FuzzyQuantizedSignal quantized;
    quantized.tokens.reserve(signal.samples.size());
    const std::int16_t sentinel = std::numeric_limits<std::int16_t>::min();

    for (const auto sample : signal.samples) {
        // Handle NaN sentinel: map to minimum int16 value
        if (std::isnan(sample)) {
            quantized.tokens.push_back(sentinel);
            continue;
        }

        const auto val = dynamic_quantize(sample, config_);
        quantized.tokens.push_back(static_cast<std::int16_t>(val));
    }
    return quantized;
}

const FuzzyQuantizerConfig& Rh2FuzzyQuantizer::config() const { return config_; }

std::string Rh2FuzzyQuantizer::name() const { return config_.backend; }

}  // namespace piru::signal
