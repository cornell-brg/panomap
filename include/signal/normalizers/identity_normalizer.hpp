// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "signal/normalizers/normalizer.hpp"

namespace piru::signal {

class IdentityNormalizer : public SignalNormalizer {
public:
    explicit IdentityNormalizer(SignalNormalizerConfig config);

    NormalizedSignal normalize(const io::RawRead& read, const EventSeries* events) const override;
    const SignalNormalizerConfig& config() const override;
    std::string name() const override;

private:
    SignalNormalizerConfig config_;
};

}  // namespace piru::signal
