// SPDX-License-Identifier: MIT

#include "signal/fuzzy_quantizers/passthrough_fuzzy_quantizer.hpp"

namespace piru::signal {

PassthroughFuzzyQuantizer::PassthroughFuzzyQuantizer(FuzzyQuantizerConfig config)
    : config_(std::move(config)) {}

FuzzyQuantizedSignal PassthroughFuzzyQuantizer::quantize(const NormalizedSignal& signal,
                                                        const EventSeries* events) const {
    (void)events;
    FuzzyQuantizedSignal quantized;
    quantized.tokens.reserve(signal.samples.size());
    for (const auto sample : signal.samples) {
        quantized.tokens.push_back(static_cast<std::int16_t>(sample));
    }
    return quantized;
}

const FuzzyQuantizerConfig& PassthroughFuzzyQuantizer::config() const { return config_; }

std::string PassthroughFuzzyQuantizer::name() const { return config_.backend; }

}  // namespace piru::signal
