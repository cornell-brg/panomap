// SPDX-License-Identifier: MIT

#include "signal/tokenizers/rh2_tokenizer.hpp"

#include <algorithm>
#include <cmath>

namespace panomap::signal {

namespace {

std::uint32_t dynamic_quantize(float value, const TokenizerConfig& cfg) {
  const float fine_min = cfg.fine_min;
  const float fine_max = cfg.fine_max;
  const float fine_range = cfg.fine_range;
  const std::uint32_t n_buckets = 1u << cfg.qbits;

  if (n_buckets == 0u) return 0u;

  float min_val = -3.0f;
  float max_val = 3.0f;
  float range = max_val - min_val;
  float coarse_coef1 = (1.0f - fine_range) * 0.5f;
  float coarse_coef2 = fine_range + coarse_coef1;

  float signal = std::clamp(value, min_val, max_val);

  // Min-max normalization to [0, 1]
  float normalized = (signal - min_val) / range;

  // Fine range mapping
  float a = (fine_min - min_val) / range;
  float b = (fine_max - min_val) / range;

  // Conditional quantization
  float quantized = fine_max;
  if (signal >= fine_min && signal <= fine_max) {
    quantized = fine_range * ((normalized - a) / (b - a));
  } else {
    if (normalized < 0.5f) {
      quantized = fine_range + coarse_coef1 * normalized;
    } else {
      quantized = coarse_coef2 + coarse_coef1 * normalized;
    }
  }

  std::uint32_t qnt_val = static_cast<std::uint32_t>(quantized * (n_buckets - 1));
  return qnt_val;
}

}  // namespace

Rh2Tokenizer::Rh2Tokenizer(TokenizerConfig config) : config_(std::move(config)) {}

TokenizedSignal Rh2Tokenizer::quantize(const NormalizedSignal& signal) const {
  TokenizedSignal result;
  result.tokens.reserve(signal.samples.size());
  const std::int16_t sentinel = std::numeric_limits<std::int16_t>::min();

  for (const auto sample : signal.samples) {
    if (std::isnan(sample)) {
      result.tokens.push_back(sentinel);
      continue;
    }
    result.tokens.push_back(static_cast<std::int16_t>(dynamic_quantize(sample, config_)));
  }

  // Pass through position mapping from diff filter (if present)
  result.original_positions = signal.original_positions;
  return result;
}

const TokenizerConfig& Rh2Tokenizer::config() const { return config_; }

std::string Rh2Tokenizer::name() const { return config_.backend; }

}  // namespace panomap::signal
