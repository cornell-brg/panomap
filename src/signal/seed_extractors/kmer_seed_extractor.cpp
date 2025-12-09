// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/kmer_seed_extractor.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace piru::signal {

namespace {

// Same mix function used in legacy pipeline; mask keeps the value bounded.
std::uint64_t hash64(std::uint64_t key, std::uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

std::uint64_t make_mask(std::uint32_t bits) {
    if (bits == 0) return std::numeric_limits<std::uint64_t>::max();
    if (bits >= 64) return std::numeric_limits<std::uint64_t>::max();
    return (1ULL << bits) - 1ULL;
}

}  // namespace

KmerSeedExtractor::KmerSeedExtractor(SeedExtractorConfig config) : config_(std::move(config)) {}

SeedBuffer KmerSeedExtractor::extract(const FuzzyQuantizedSignal& signal,
                                      const EventSeries* events) const {
    (void)events;
    SeedBuffer buffer;
    const auto& tokens = signal.tokens;
    const std::size_t k = config_.k;
    const std::size_t stride = std::max<std::size_t>(config_.stride, 1);

    if (k == 0 || tokens.size() < k) return buffer;

    const std::uint32_t qbits = config_.qbits;
    const std::uint64_t token_mask = make_mask(qbits);

    // Guard against oversized windows (shift by >=64 is UB). Fall back to full mask.
    const bool use_shift = qbits > 0 && qbits < 64;
    const std::uint64_t window_mask =
        (!use_shift || k >= (64 / static_cast<std::size_t>(qbits)))
            ? std::numeric_limits<std::uint64_t>::max()
            : make_mask(static_cast<std::uint32_t>(k * qbits));

    const std::size_t num_seeds = 1 + (tokens.size() - k) / stride;
    buffer.seeds.reserve(num_seeds);

    // Note: k-mer backend produces 32-bit hashes (legacy compatibility) stored in 64-bit
    // field. The 64-bit Seed.hash field allows future backends (e.g., minimizer) to use
    // full 64-bit hash space if needed.
    const std::uint64_t mix_mask = 0xFFFFFFFFULL;

    const std::int16_t sentinel = std::numeric_limits<std::int16_t>::min();

    for (std::size_t start = 0; start + k <= tokens.size(); start += stride) {
        // Check if window contains sentinel (N base marker)
        bool has_sentinel = false;
        for (std::size_t j = 0; j < k; ++j) {
            if (tokens[start + j] == sentinel) {
                has_sentinel = true;
                break;
            }
        }

        // Skip seed extraction for windows containing sentinels
        if (has_sentinel) {
            continue;
        }

        std::uint64_t packed = 0;
        // Two packing strategies to avoid undefined behavior:
        // 1. Normal case (k*qbits < 64): bit-pack tokens via shift-and-OR.
        // 2. Overflow case (k*qbits >= 64 or qbits==0): incrementally mix tokens via
        //    hash function to prevent shift-by-64+ undefined behavior.
        for (std::size_t j = 0; j < k; ++j) {
            const auto token = static_cast<std::uint64_t>(
                static_cast<std::uint16_t>(tokens[start + j])) &
                              token_mask;
            packed = (!use_shift)
                         ? hash64(packed ^ token, mix_mask)  // Overflow: mix via hash
                         : ((packed << qbits) | token) & window_mask;  // Normal: pack bits
        }
        const auto hash = hash64(packed, mix_mask);
        buffer.seeds.push_back(Seed{.hash = hash, .position = start, .span = k});
    }
    return buffer;
}

const SeedExtractorConfig& KmerSeedExtractor::config() const { return config_; }

std::string KmerSeedExtractor::name() const { return config_.backend; }

}  // namespace piru::signal
