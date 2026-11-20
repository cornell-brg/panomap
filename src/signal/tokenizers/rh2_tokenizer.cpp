// SPDX-License-Identifier: MIT

#include "signal/tokenizers/rh2_tokenizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace panomap::signal {

namespace {

// N ordered buckets: the fine range gets round(N * fine_range) of them, the tails
// below and above share the rest. Matches RawHash2 2.5 dynamic_quantize.
std::uint32_t bucket_quantize(float value, const TokenizerConfig& cfg) {
  const std::uint32_t n = cfg.n_buckets;
  const std::uint32_t n_fine =
      std::clamp<std::uint32_t>(static_cast<std::uint32_t>(n * cfg.fine_range + 0.5f), 1u, n);
  const std::uint32_t n_low = (n - n_fine) / 2;
  const std::uint32_t n_high = n - n_fine - n_low;
  const float lo = -3.0f;
  const float hi = 3.0f;
  auto bucket = [](float t, std::uint32_t count) {
    t = std::clamp(t, 0.0f, 1.0f);
    return std::min(static_cast<std::uint32_t>(t * count), count - 1);
  };
  if (value < cfg.fine_min) {
    if (n_low == 0) return 0;
    return bucket((value - lo) / (cfg.fine_min - lo), n_low);
  }
  if (value <= cfg.fine_max) {
    return n_low + bucket((value - cfg.fine_min) / (cfg.fine_max - cfg.fine_min), n_fine);
  }
  if (n_high == 0) return n - 1;
  return n_low + n_fine + bucket((value - cfg.fine_max) / (hi - cfg.fine_max), n_high);
}

std::uint32_t dynamic_quantize(float value, const TokenizerConfig& cfg) {
  if (cfg.n_buckets > 0) return bucket_quantize(value, cfg);
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

// Drop each token equal to the last kept token. A sentinel (NaN sample) is always
// kept and resets the comparison.
// Kept tokens keep their original event positions, so seed positions stay in
// event units.
void collapse_repeats(TokenizedSignal& signal, std::int16_t sentinel) {
  if (signal.tokens.empty()) return;
  const bool has_pos = !signal.original_positions.empty();
  std::vector<std::int16_t> kept;
  std::vector<std::uint32_t> kept_pos;
  kept.reserve(signal.tokens.size());
  kept_pos.reserve(signal.tokens.size());
  bool have_prev = false;
  std::int16_t prev = 0;
  for (std::size_t i = 0; i < signal.tokens.size(); ++i) {
    const std::int16_t tok = signal.tokens[i];
    if (tok != sentinel && have_prev && tok == prev) continue;
    kept.push_back(tok);
    kept_pos.push_back(has_pos ? signal.original_positions[i] : static_cast<std::uint32_t>(i));
    have_prev = tok != sentinel;
    prev = tok;
  }
  signal.tokens = std::move(kept);
  signal.original_positions = std::move(kept_pos);
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
  if (config_.collapse_repeats) collapse_repeats(result, sentinel);
  return result;
}

const TokenizerConfig& Rh2Tokenizer::config() const { return config_; }

std::string Rh2Tokenizer::name() const { return config_.backend; }

}  // namespace panomap::signal
