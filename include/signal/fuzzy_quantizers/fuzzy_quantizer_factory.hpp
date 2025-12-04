// SPDX-License-Identifier: MIT
// Factory helpers for fuzzy quantizers.

#pragma once

#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"

namespace piru::signal {

FuzzyQuantizerPtr make_fuzzy_quantizer(const FuzzyQuantizerConfig& config);

}  // namespace piru::signal
