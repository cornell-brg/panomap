// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/minimizer_seed_extractor.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "signal/seed_extractors/seed_hash_util.hpp"

namespace piru::signal {

MinimizerSeedExtractor::MinimizerSeedExtractor(SeedExtractorConfig config)
    : config_(std::move(config)) {}

SeedBuffer MinimizerSeedExtractor::extract(const TokenizedSignal& signal) const {
  SeedBuffer buffer;
  const auto& tokens = signal.tokens;
  const std::size_t k = config_.k;
  const std::size_t w = std::max<std::size_t>(config_.window, 1);

  if (k == 0 || tokens.size() < k) return buffer;

  const std::uint32_t qbits = config_.qbits;
  const std::uint64_t token_mask = makeMask(qbits);

  const bool use_shift = qbits > 0 && qbits < 64;
  const std::uint64_t window_mask = (!use_shift || k >= (64 / static_cast<std::size_t>(qbits)))
                                        ? std::numeric_limits<std::uint64_t>::max()
                                        : makeMask(static_cast<std::uint32_t>(k * qbits));

  // Step 1: Compute all k-mer hashes.
  const std::size_t num_kmers = tokens.size() - k + 1;
  constexpr std::uint64_t kInvalidHash = std::numeric_limits<std::uint64_t>::max();
  std::vector<std::uint64_t> hashes(num_kmers, kInvalidHash);

  for (std::size_t i = 0; i < num_kmers; ++i) {
    if (hasSentinel(tokens.data(), i, k)) continue;
    hashes[i] = hashKmer(tokens.data(), i, k, qbits, token_mask, window_mask, use_shift);
  }

  // Step 2: Slide window of size w, emit minimizer (leftmost on tie).
  if (num_kmers < w) return buffer;

  buffer.seeds.reserve(num_kmers / w + 1);
  std::size_t last_emitted = std::numeric_limits<std::size_t>::max();

  for (std::size_t i = 0; i <= num_kmers - w; ++i) {
    // Find leftmost minimum in window [i, i+w).
    std::size_t min_idx = i;
    std::uint64_t min_hash = hashes[i];
    for (std::size_t j = i + 1; j < i + w; ++j) {
      if (hashes[j] < min_hash) {
        min_hash = hashes[j];
        min_idx = j;
      }
    }

    // Skip windows where the minimizer is a sentinel.
    if (min_hash == kInvalidHash) continue;

    // Deduplicate: don't emit the same position twice.
    if (min_idx != last_emitted) {
      const bool has_pos = !signal.original_positions.empty();
      const std::size_t pos = has_pos ? signal.original_positions[min_idx] : min_idx;
      const std::size_t end = has_pos ? signal.original_positions[min_idx + k - 1] : min_idx + k - 1;
      const std::size_t span = end - pos + 1;
      buffer.seeds.push_back(Seed{.hash = min_hash, .position = pos, .length = span});
      last_emitted = min_idx;
    }
  }

  return buffer;
}

const SeedExtractorConfig& MinimizerSeedExtractor::config() const { return config_; }

std::string MinimizerSeedExtractor::name() const { return "minimizer"; }

}  // namespace piru::signal
