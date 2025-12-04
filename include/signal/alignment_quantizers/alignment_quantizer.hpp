// SPDX-License-Identifier: MIT
// Interface for quantizing signals into alignment-ready representations.

#pragma once

#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct AlignmentQuantizerConfig {
    std::string backend{"int16"};
    std::string params;
    bool operator==(const AlignmentQuantizerConfig& other) const = default;
};

class AlignmentQuantizer {
public:
    virtual ~AlignmentQuantizer() = default;

    virtual AlignmentQuantizedSignal quantize(const NormalizedSignal& signal,
                                              const EventSeries* events) const = 0;
    virtual const AlignmentQuantizerConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using AlignmentQuantizerPtr = std::unique_ptr<AlignmentQuantizer>;

}  // namespace piru::signal
