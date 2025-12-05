// SPDX-License-Identifier: MIT

#include "signal/normalizers/identity_normalizer.hpp"

namespace piru::signal {

IdentityNormalizer::IdentityNormalizer(SignalNormalizerConfig config)
    : config_(std::move(config)) {}

NormalizedSignal IdentityNormalizer::normalize(const io::RawRead& read, const EventSeries* events) const {
    (void)events;  // Unused for now
    NormalizedSignal normalized;
    normalized.samples.reserve(read.raw_signal.size());
    for (const auto value : read.raw_signal) {
        normalized.samples.push_back(static_cast<float>(value));
    }
    normalized.sampling_rate_hz = read.sampling_rate_hz;
    return normalized;
}

const SignalNormalizerConfig& IdentityNormalizer::config() const { return config_; }

std::string IdentityNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
