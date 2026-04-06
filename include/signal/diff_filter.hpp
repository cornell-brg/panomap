// SPDX-License-Identifier: MIT
// Diff filter: remove near-duplicate events from a normalized signal.

#pragma once

#include "signal/signal_types.hpp"

namespace piru::signal {

// Drop events whose value is within `diff` of the last emitted event.
// Compresses the signal in place and populates original_positions to map
// surviving indices back to the original event stream.
// No-op when diff <= 0.
void apply_diff_filter(NormalizedSignal& signal, float diff);

}  // namespace piru::signal
