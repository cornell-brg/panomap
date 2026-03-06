// SPDX-License-Identifier: MIT

#include "signal/fuzzy_quantizers/piru_fuzzy_quantizer.hpp"

#include <algorithm>
#include <cmath>

namespace piru::signal {

namespace {

/// Quantize with proper coarse-fine-coarse distribution using ALL buckets.
///
/// Layout:
///   [min_val, fine_min) -> first (1-fine_range)/2 of buckets  (lower coarse)
///   [fine_min, fine_max] -> middle fine_range of buckets      (fine)
///   (fine_max, max_val] -> last (1-fine_range)/2 of buckets   (upper coarse)
std::uint32_t piru_quantize(float value, const FuzzyQuantizerConfig& cfg) {
  const float fine_min = cfg.fine_min;
  const float fine_max = cfg.fine_max;
  const float fine_range = cfg.fine_range;
  const std::uint32_t n_buckets = cfg.n_bins > 0 ? cfg.n_bins : (1u << cfg.qbits);

  if (n_buckets == 0u) return 0u;

  constexpr float min_val = -3.0f;
  constexpr float max_val = 3.0f;

  // Clamp signal to valid range
  float signal = std::clamp(value, min_val, max_val);

  // Calculate bucket boundaries
  const float coarse_range = (1.0f - fine_range) / 2.0f;
  const float lower_coarse_end = coarse_range;       // e.g., 0.3
  const float fine_end = coarse_range + fine_range;  // e.g., 0.7
  // upper_coarse_end = 1.0

  float quantized;

  if (signal < fine_min) {
    // Lower coarse region: [min_val, fine_min) -> [0, coarse_range)
    float t = (signal - min_val) / (fine_min - min_val);  // [0, 1)
    quantized = t * coarse_range;
  } else if (signal <= fine_max) {
    // Fine region: [fine_min, fine_max] -> [coarse_range, coarse_range + fine_range]
    float t = (signal - fine_min) / (fine_max - fine_min);  // [0, 1]
    quantized = lower_coarse_end + t * fine_range;
  } else {
    // Upper coarse region: (fine_max, max_val] -> (fine_end, 1.0]
    float t = (signal - fine_max) / (max_val - fine_max);  // (0, 1]
    quantized = fine_end + t * coarse_range;
  }

  // Map [0, 1] to [0, n_buckets-1]
  std::uint32_t bucket = static_cast<std::uint32_t>(quantized * (n_buckets - 1));
  return std::min(bucket, n_buckets - 1);
}

}  // namespace

PiruFuzzyQuantizer::PiruFuzzyQuantizer(FuzzyQuantizerConfig config) : config_(std::move(config)) {}

FuzzyQuantizedSignal PiruFuzzyQuantizer::quantize(const NormalizedSignal& signal) const {
  FuzzyQuantizedSignal quantized;
  quantized.tokens.reserve(signal.samples.size());
  const std::int16_t sentinel = std::numeric_limits<std::int16_t>::min();

  for (const auto sample : signal.samples) {
    // Handle NaN sentinel: map to minimum int16 value
    if (std::isnan(sample)) {
      quantized.tokens.push_back(sentinel);
      continue;
    }

    const auto val = piru_quantize(sample, config_);
    quantized.tokens.push_back(static_cast<std::int16_t>(val));
  }
  return quantized;
}

const FuzzyQuantizerConfig& PiruFuzzyQuantizer::config() const { return config_; }

std::string PiruFuzzyQuantizer::name() const { return "piru"; }

}  // namespace piru::signal
