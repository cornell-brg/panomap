// SPDX-License-Identifier: MIT

#include "signal/tokenizers/tokenizer_factory.hpp"

#include "signal/tokenizers/landmark_tokenizer.hpp"
#include "signal/tokenizers/rh2_tokenizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

TokenizerPtr make_tokenizer(const TokenizerConfig& config) {
  if (config.backend == "rawhash2" || config.backend == "rh2" || config.backend.empty()) {
    return std::make_unique<Rh2Tokenizer>(config);
  }
  if (config.backend == "landmark") {
    LandmarkTokenizerConfig lm_cfg;
    lm_cfg.rise_bits = 2;
    lm_cfg.fall_bits = 2;
    lm_cfg.max_amp = 4.0f;
    lm_cfg.min_prominence = config.landmark_min_prominence;
    return std::make_unique<LandmarkTokenizer>(config, lm_cfg);
  }

  LOG_WARN("Unknown tokenizer backend '" + config.backend + "', using rh2");
  return std::make_unique<Rh2Tokenizer>(config);
}

}  // namespace piru::signal
