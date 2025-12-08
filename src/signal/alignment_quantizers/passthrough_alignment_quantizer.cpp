// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/passthrough_alignment_quantizer.hpp"

#include <vector>

namespace piru::signal {

PassthroughAlignmentQuantizer::PassthroughAlignmentQuantizer(AlignmentQuantizerConfig config)
    : config_(std::move(config)) {}

AlignmentQuantizedSignal PassthroughAlignmentQuantizer::quantize(const NormalizedSignal& signal,
                                                                 const EventSeries* events) const {
    (void)events;
    AlignmentQuantizedSignal quantized;
    quantized.kind = AlignmentQuantizationKind::kInt16;
    std::vector<std::int16_t> values;
    values.reserve(signal.samples.size());
    for (const auto sample : signal.samples) {
        values.push_back(static_cast<std::int16_t>(sample));
    }
    quantized.data = std::move(values);
    return quantized;
}

const AlignmentQuantizerConfig& PassthroughAlignmentQuantizer::config() const {
    return config_;
}

std::string PassthroughAlignmentQuantizer::name() const { return config_.backend; }

float PassthroughAlignmentQuantizer::scale() const { return 1.0f; }

float PassthroughAlignmentQuantizer::offset() const { return 0.0f; }

}  // namespace piru::signal
