// SPDX-License-Identifier: MIT

#pragma once

#include "signal/normalizers/normalizer.hpp"

namespace piru::signal {

class ZScoreNormalizer : public SignalNormalizer {
public:
    explicit ZScoreNormalizer(SignalNormalizerConfig config);

    NormalizedSignal normalize(const io::RawRead& read, const EventSeries* events) const override;
    const SignalNormalizerConfig& config() const override;
    std::string name() const override;

private:
    SignalNormalizerConfig config_;
};

}  // namespace piru::signal
