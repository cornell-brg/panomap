// SPDX-License-Identifier: MIT

#pragma once

#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"

namespace piru::signal {

class Rh2FuzzyQuantizer : public FuzzyQuantizer {
public:
    explicit Rh2FuzzyQuantizer(FuzzyQuantizerConfig config);

    FuzzyQuantizedSignal quantize(const NormalizedSignal& signal,
                                  const EventSeries* events) const override;
    const FuzzyQuantizerConfig& config() const override;
    std::string name() const override;

private:
    FuzzyQuantizerConfig config_;
};

}  // namespace piru::signal
