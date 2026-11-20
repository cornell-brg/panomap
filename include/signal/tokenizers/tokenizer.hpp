// SPDX-License-Identifier: MIT
// Interface for signal tokenization used during seeding.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace panomap::signal {

struct TokenizerConfig {
  std::string backend{"rh2"};
  std::string params;
  std::string pore_model;  // For chemistry-specific defaults (r9/r10)
  float fine_min{-2.0f};
  float fine_max{2.0f};
  float fine_range{0.85f};  // Original default: 0.4f (DEV025: R9=0.9, R10=0.8)
  std::uint32_t qbits{4};
  // 0: the original RawHash2 formula over 2^qbits values. N > 0: N ordered buckets,
  // with round(N * fine_range) of them in [fine_min, fine_max] and the rest split
  // between the two tails (RawHash2 2.5).
  std::uint32_t n_buckets{0};
  float landmark_min_prominence{0.5f};  // Landmark: drop peaks with prominence below this
  std::uint32_t landmark_bits{2};       // Landmark: bits for each of rise and fall
  float landmark_max_amp{4.0f};         // Landmark: amplitude clip for log quantization
  // Drop a token equal to the last kept token. Kept tokens keep their original
  // event positions. Recorded in the index.
  bool collapse_repeats{false};
  bool operator==(const TokenizerConfig& other) const = default;
};

// Bits per token: 2 * landmark_bits for landmark, ceil(log2(n_buckets)) for rh2
// buckets, else 4.
inline std::uint32_t token_bits(const TokenizerConfig& cfg) {
  if (cfg.backend == "landmark") return 2 * cfg.landmark_bits;
  if (cfg.n_buckets > 16) {
    std::uint32_t q = 0;
    while ((1u << q) < cfg.n_buckets) ++q;
    return q;
  }
  return 4;
}

class Tokenizer {
public:
  virtual ~Tokenizer() = default;

  virtual TokenizedSignal quantize(const NormalizedSignal& signal) const = 0;
  virtual const TokenizerConfig& config() const = 0;
  virtual std::string name() const = 0;
};

using TokenizerPtr = std::unique_ptr<Tokenizer>;

}  // namespace panomap::signal
