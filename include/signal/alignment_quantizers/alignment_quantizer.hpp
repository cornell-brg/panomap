// SPDX-License-Identifier: MIT
// Interface for quantizing signals into alignment-ready representations.

#pragma once

#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct AlignmentQuantizerConfig {
    std::string backend{"passthrough"};
    int target_bits{16};
    float scale{1.0f};
    float offset{0.0f};
};

class AlignmentQuantizer {

public:

    virtual ~AlignmentQuantizer() = default;



    virtual AlignmentQuantizedSignal quantize(const NormalizedSignal& signal,

                                              const EventSeries* events) const = 0;



    virtual const AlignmentQuantizerConfig& config() const = 0;

    virtual std::string name() const = 0;

    virtual float scale() const = 0;

    virtual float offset() const = 0;

};

using AlignmentQuantizerPtr = std::unique_ptr<AlignmentQuantizer>;

}  // namespace piru::signal
