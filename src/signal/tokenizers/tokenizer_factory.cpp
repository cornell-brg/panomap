// SPDX-License-Identifier: MIT

#include "signal/tokenizers/tokenizer_factory.hpp"

#include "signal/tokenizers/rh2_tokenizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

TokenizerPtr make_tokenizer(const TokenizerConfig& config) {
  if (config.backend == "rawhash2" || config.backend == "rh2" || config.backend.empty()) {
    return std::make_unique<Rh2Tokenizer>(config);
  }

  LOG_WARN("Unknown tokenizer backend '" + config.backend + "', using rh2");
  return std::make_unique<Rh2Tokenizer>(config);
}

}  // namespace piru::signal
