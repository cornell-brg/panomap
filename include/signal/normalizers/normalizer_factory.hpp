// SPDX-License-Identifier: MIT
// Factory helpers for signal normalizers.

#pragma once

#include "signal/normalizers/normalizer.hpp"

namespace piru::signal {

SignalNormalizerPtr make_signal_normalizer(const SignalNormalizerConfig& config);

}  // namespace piru::signal
