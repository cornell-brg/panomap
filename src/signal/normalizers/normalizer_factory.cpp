// SPDX-License-Identifier: MIT

#include "signal/normalizers/normalizer_factory.hpp"

#include "signal/normalizers/identity_normalizer.hpp"
#include "signal/normalizers/median_mad_normalizer.hpp"
#include "signal/normalizers/zscore_normalizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

SignalNormalizerPtr make_signal_normalizer(const SignalNormalizerConfig& config) {
    if (config.backend == "zscore" || config.backend.empty()) {
        return std::make_unique<ZScoreNormalizer>(config);
    }
    if (config.backend == "median_mad") {
        return std::make_unique<MedianMadNormalizer>(config);
    }
    if (config.backend != "identity") {
        LOG_WARN("Unknown signal normalizer backend '" + config.backend + "', using identity");
    }
    return std::make_unique<IdentityNormalizer>(config);
}

}  // namespace piru::signal
