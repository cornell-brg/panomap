// SPDX-License-Identifier: MIT
// Factory helpers for tokenizers.

#pragma once

#include "signal/tokenizers/tokenizer.hpp"

namespace piru::signal {

TokenizerPtr make_tokenizer(const TokenizerConfig& config);

}  // namespace piru::signal
