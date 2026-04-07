// SPDX-License-Identifier: MIT
// Landmark tokenizer: encode signal peaks by rise/fall amplitudes.
//
// Detects valleys (landmarks) in the normalized signal, extracts peaks
// between consecutive valleys, and encodes each peak as a token from
// its log-quantized rise and fall amplitudes.
//
// Related:
//  - python-rh2/pyrh2/peak_fingerprint.py (prototype)
//  - python-rh2/docs/peakseed.md (design doc)

#pragma once

#include "signal/tokenizers/tokenizer.hpp"

namespace piru::signal {

struct LandmarkTokenizerConfig {
  std::uint32_t rise_bits{2};      // bits for rise amplitude (4 levels)
  std::uint32_t fall_bits{2};      // bits for fall amplitude (4 levels)
  float max_amp{4.0f};             // clip amplitude for log quantization
  float min_prominence{0.5f};      // minimum peak prominence to keep
};

class LandmarkTokenizer : public Tokenizer {
public:
  explicit LandmarkTokenizer(TokenizerConfig config, LandmarkTokenizerConfig landmark_config);

  TokenizedSignal quantize(const NormalizedSignal& signal) const override;
  const TokenizerConfig& config() const override;
  std::string name() const override;

private:
  TokenizerConfig config_;
  LandmarkTokenizerConfig landmark_config_;
};

}  // namespace piru::signal
