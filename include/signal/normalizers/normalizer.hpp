// SPDX-License-Identifier: MIT
// Interface for normalizing raw read signals.

#pragma once

#include <memory>
#include <string>

#include "io/reads/read_provider.hpp"
#include "signal/signal_types.hpp"

namespace piru::signal {

struct SignalNormalizerConfig {
    std::string backend{"zscore"};
};

class SignalNormalizer {
public:
    virtual ~SignalNormalizer() = default;

    virtual NormalizedSignal normalize(const io::RawRead& read) const = 0;
    virtual const SignalNormalizerConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using SignalNormalizerPtr = std::unique_ptr<SignalNormalizer>;

}  // namespace piru::signal
