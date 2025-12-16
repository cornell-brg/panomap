// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/passthrough_alignment_quantizer.hpp"

#include <vector>

namespace piru::signal {

PassthroughAlignmentQuantizer::PassthroughAlignmentQuantizer(AlignmentQuantizerConfig config)
    : config_(std::move(config)) {}

AlignmentQuantizedSignal PassthroughAlignmentQuantizer::quantize(const NormalizedSignal& signal) const {
    AlignmentQuantizedSignal quantized;
    quantized.kind = AlignmentQuantizationKind::kFloat32;
    quantized.data = signal.samples;
    return quantized;
}

const AlignmentQuantizerConfig& PassthroughAlignmentQuantizer::config() const {
    return config_;
}

std::string PassthroughAlignmentQuantizer::name() const { return config_.backend; }

float PassthroughAlignmentQuantizer::scale() const { return 1.0f; }

float PassthroughAlignmentQuantizer::offset() const { return 0.0f; }

}  // namespace piru::signal
