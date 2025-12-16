// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/int_alignment_quantizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace piru::signal {

namespace {

template <typename T>
std::vector<T> quantize_to_int(const std::vector<float>& samples, float scale) {
    std::vector<T> out;
    out.reserve(samples.size());
    const auto min_val = static_cast<double>(std::numeric_limits<T>::min());
    const auto max_val = static_cast<double>(std::numeric_limits<T>::max());
    const T sentinel = std::numeric_limits<T>::min();

    for (const auto s : samples) {
        // Handle NaN sentinel: map to minimum value
        if (std::isnan(s)) {
            out.push_back(sentinel);
            continue;
        }

        double v = static_cast<double>(s) * static_cast<double>(scale);
        v = std::clamp(std::round(v), min_val + 1.0, max_val);  // Reserve min for sentinel
        out.push_back(static_cast<T>(v));
    }
    return out;
}

}  // namespace

IntAlignmentQuantizer::IntAlignmentQuantizer(AlignmentQuantizerConfig config)
    : config_(std::move(config)) {
    // If the scale is the default, calculate a new one to maximize the dynamic range
    // based on the normalized input range of [-3.0, 3.0].
    if (config_.scale == 1.0f) {
        const float input_max = 3.0f;
        if (config_.target_bits == 8) {
            // Reserve min for sentinel
            config_.scale = (static_cast<float>(std::numeric_limits<int8_t>::max()) - 1.0f) / input_max;
        } else {  // Default to 16 bits
            config_.target_bits = 16;
            // Reserve min for sentinel
            config_.scale = (static_cast<float>(std::numeric_limits<int16_t>::max()) - 1.0f) / input_max;
        }
    }
}

AlignmentQuantizedSignal IntAlignmentQuantizer::quantize(const NormalizedSignal& signal) const {
    AlignmentQuantizedSignal quantized;
    const int bits = (config_.target_bits == 8) ? 8 : 16;
    
    if (bits == 8) {
        quantized.kind = AlignmentQuantizationKind::kInt8;
        quantized.data = quantize_to_int<std::int8_t>(signal.samples, config_.scale);
    } else {
        quantized.kind = AlignmentQuantizationKind::kInt16;
        quantized.data = quantize_to_int<std::int16_t>(signal.samples, config_.scale);
    }
    return quantized;
}

const AlignmentQuantizerConfig& IntAlignmentQuantizer::config() const { return config_; }

std::string IntAlignmentQuantizer::name() const { return config_.backend; }

float IntAlignmentQuantizer::scale() const { return config_.scale; }

float IntAlignmentQuantizer::offset() const { return config_.offset; }

}  // namespace piru::signal
