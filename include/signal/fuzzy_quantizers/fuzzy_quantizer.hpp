// SPDX-License-Identifier: MIT
// Interface for fuzzy quantization used during seeding.

#pragma once

#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct FuzzyQuantizerConfig {
    std::string backend{"rawhash2"};
    std::string params;
    float fine_min{-2.0f};
    float fine_max{2.0f};
    float fine_range{0.4f};
    std::uint32_t qbits{4};
    bool operator==(const FuzzyQuantizerConfig& other) const = default;
};

class FuzzyQuantizer {
public:
    virtual ~FuzzyQuantizer() = default;

    virtual FuzzyQuantizedSignal quantize(const NormalizedSignal& signal,
                                          const EventSeries* events) const = 0;
    virtual const FuzzyQuantizerConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using FuzzyQuantizerPtr = std::unique_ptr<FuzzyQuantizer>;

}  // namespace piru::signal
