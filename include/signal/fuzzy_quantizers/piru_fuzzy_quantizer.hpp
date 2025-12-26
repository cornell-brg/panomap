// SPDX-License-Identifier: MIT

#pragma once

#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"

namespace piru::signal {

/// Piru fuzzy quantizer with proper coarse-fine-coarse bucket distribution.
/// Unlike RH2, this uses ALL available buckets with no gaps.
///
/// Bucket layout (e.g., 16 bins, fine_range=0.4):
///   Bins 0-4   (30%): Lower coarse [min_val, fine_min)
///   Bins 5-10  (40%): Fine region  [fine_min, fine_max]
///   Bins 11-15 (30%): Upper coarse (fine_max, max_val]
class PiruFuzzyQuantizer : public FuzzyQuantizer {
public:
    explicit PiruFuzzyQuantizer(FuzzyQuantizerConfig config);

    FuzzyQuantizedSignal quantize(const NormalizedSignal& signal) const override;
    const FuzzyQuantizerConfig& config() const override;
    std::string name() const override;

private:
    FuzzyQuantizerConfig config_;
};

}  // namespace piru::signal
