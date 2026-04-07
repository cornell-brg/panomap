// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/kmer_seed_extractor.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

#include "signal/seed_extractors/seed_hash_util.hpp"

namespace piru::signal {

KmerSeedExtractor::KmerSeedExtractor(SeedExtractorConfig config) : config_(std::move(config)) {}

SeedBuffer KmerSeedExtractor::extract(const TokenizedSignal& signal) const {
  SeedBuffer buffer;
  const auto& tokens = signal.tokens;
  const std::size_t k = config_.k;
  const std::size_t stride = std::max<std::size_t>(config_.stride, 1);

  if (k == 0 || tokens.size() < k) return buffer;

  const std::uint32_t qbits = config_.qbits;
  const std::uint64_t token_mask = makeMask(qbits);

  // Guard against oversized windows (shift by >=64 is UB). Fall back to full mask.
  const bool use_shift = qbits > 0 && qbits < 64;
  const std::uint64_t window_mask = (!use_shift || k >= (64 / static_cast<std::size_t>(qbits)))
                                        ? std::numeric_limits<std::uint64_t>::max()
                                        : makeMask(static_cast<std::uint32_t>(k * qbits));

  const std::size_t num_seeds = 1 + (tokens.size() - k) / stride;
  buffer.seeds.reserve(num_seeds);

  for (std::size_t start = 0; start + k <= tokens.size(); start += stride) {
    if (hasSentinel(tokens.data(), start, k)) continue;

    const auto hash = hashKmer(tokens.data(), start, k, qbits, token_mask, window_mask, use_shift);
    const bool has_pos = !signal.original_positions.empty();
    const std::size_t pos = has_pos ? signal.original_positions[start] : start;
    const std::size_t end = has_pos ? signal.original_positions[start + k - 1] : start + k - 1;
    const std::size_t span = end - pos + 1;
    buffer.seeds.push_back(Seed{.hash = hash, .position = pos, .length = span});
  }
  return buffer;
}

const SeedExtractorConfig& KmerSeedExtractor::config() const { return config_; }

std::string KmerSeedExtractor::name() const { return config_.backend; }

}  // namespace piru::signal
