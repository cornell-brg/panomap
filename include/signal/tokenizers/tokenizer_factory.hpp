// SPDX-License-Identifier: MIT
// Factory helpers for tokenizers.

#pragma once

#include "signal/tokenizers/tokenizer.hpp"

namespace panomap::signal {

TokenizerPtr make_tokenizer(const TokenizerConfig& config);

}  // namespace panomap::signal
