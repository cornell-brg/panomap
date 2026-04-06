// SPDX-License-Identifier: MIT

#pragma once

#include "signal/tokenizers/tokenizer.hpp"

namespace piru::signal {

class Rh2Tokenizer : public Tokenizer {
public:
  explicit Rh2Tokenizer(TokenizerConfig config);

  TokenizedSignal quantize(const NormalizedSignal& signal) const override;
  const TokenizerConfig& config() const override;
  std::string name() const override;

private:
  TokenizerConfig config_;
};

}  // namespace piru::signal
