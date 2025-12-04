// SPDX-License-Identifier: MIT
// Factory helpers for alignment quantizers.

#pragma once

#include "signal/alignment_quantizers/alignment_quantizer.hpp"

namespace piru::signal {

AlignmentQuantizerPtr make_alignment_quantizer(const AlignmentQuantizerConfig& config);

}  // namespace piru::signal
