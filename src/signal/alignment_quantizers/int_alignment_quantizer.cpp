// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/int_alignment_quantizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace piru::signal {

namespace {

template <typename T>
std::vector<T> quantize_to_int(const std::vector<float>& samples, float scale) {
    std::vector<T> out;
    out.reserve(samples.size());
    const auto min_val = static_cast<double>(std::numeric_limits<T>::min());
    const auto max_val = static_cast<double>(std::numeric_limits<T>::max());
    for (const auto s : samples) {
        double v = static_cast<double>(s) * static_cast<double>(scale);
        v = std::clamp(std::round(v), min_val, max_val);
        out.push_back(static_cast<T>(v));
    }
    return out;
}

}  // namespace

IntAlignmentQuantizer::IntAlignmentQuantizer(AlignmentQuantizerConfig config)
    : config_(std::move(config)) {}

AlignmentQuantizedSignal IntAlignmentQuantizer::quantize(const NormalizedSignal& signal,
                                                         const EventSeries* events) const {
    (void)events;
    AlignmentQuantizedSignal quantized;
    const int bits = (config_.target_bits == 8) ? 8 : 16;
    const float scale = 1.0f;  // Placeholder for future fixed-point scaling tweaks.

    if (bits == 8) {
        quantized.kind = AlignmentQuantizationKind::kInt8;
        quantized.data = quantize_to_int<std::int8_t>(signal.samples, scale);
    } else {
        quantized.kind = AlignmentQuantizationKind::kInt16;
        quantized.data = quantize_to_int<std::int16_t>(signal.samples, scale);
    }
    return quantized;
}

const AlignmentQuantizerConfig& IntAlignmentQuantizer::config() const { return config_; }

std::string IntAlignmentQuantizer::name() const { return config_.backend; }

}  // namespace piru::signal
