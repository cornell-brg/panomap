// SPDX-License-Identifier: MIT

#include "signal/normalizers/identity_normalizer.hpp"

namespace piru::signal {

IdentityNormalizer::IdentityNormalizer(SignalNormalizerConfig config)
    : config_(std::move(config)) {}

NormalizedSignal IdentityNormalizer::normalize(const EventSeries& events) const {
    NormalizedSignal normalized;
    normalized.sampling_rate_hz = events.sampling_rate_hz;

    // Pass through event means (no normalization)
    normalized.samples.reserve(events.events.size());
    for (const auto& event : events.events) {
        normalized.samples.push_back(event.mean);
    }
    return normalized;
}

const SignalNormalizerConfig& IdentityNormalizer::config() const { return config_; }

std::string IdentityNormalizer::name() const { return config_.backend; }

}  // namespace piru::signal
