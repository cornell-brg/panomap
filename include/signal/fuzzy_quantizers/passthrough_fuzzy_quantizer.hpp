// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"

namespace piru::signal {

class PassthroughFuzzyQuantizer : public FuzzyQuantizer {
public:
    explicit PassthroughFuzzyQuantizer(FuzzyQuantizerConfig config);

    FuzzyQuantizedSignal quantize(const NormalizedSignal& signal,
                                  const EventSeries* events) const override;
    const FuzzyQuantizerConfig& config() const override;
    std::string name() const override;

private:
    FuzzyQuantizerConfig config_;
};

}  // namespace piru::signal
