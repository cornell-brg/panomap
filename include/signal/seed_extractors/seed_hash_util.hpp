// SPDX-License-Identifier: MIT
// Shared hashing utilities for seed extractors.

#pragma once

#include <cstdint>
#include <limits>

namespace piru::signal {

// Mix function used across seed backends; mask keeps the value bounded.
inline std::uint64_t hash64(std::uint64_t key, std::uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

inline std::uint64_t makeMask(std::uint32_t bits) {
    if (bits == 0 || bits >= 64) return std::numeric_limits<std::uint64_t>::max();
    return (1ULL << bits) - 1ULL;
}

// Hash mix mask for 32-bit output (legacy compatibility).
inline constexpr std::uint64_t kMixMask = 0xFFFFFFFFULL;

// Sentinel value marking N-bases in quantized token streams.
inline constexpr std::int16_t kTokenSentinel = std::numeric_limits<std::int16_t>::min();

// Compute the hash for a k-mer starting at tokens[start].
// Returns the hash, or std::nullopt-equivalent behavior is handled by the caller
// via the has_sentinel check.
inline std::uint64_t hashKmer(const std::int16_t* tokens, std::size_t start,
                              std::size_t k, std::uint32_t qbits,
                              std::uint64_t token_mask, std::uint64_t window_mask,
                              bool use_shift) {
    std::uint64_t packed = 0;
    for (std::size_t j = 0; j < k; ++j) {
        const auto token =
            static_cast<std::uint64_t>(static_cast<std::uint16_t>(tokens[start + j])) & token_mask;
        packed = (!use_shift) ? hash64(packed ^ token, kMixMask)
                              : ((packed << qbits) | token) & window_mask;
    }
    return hash64(packed, kMixMask);
}

// Check if a k-mer window starting at tokens[start] contains a sentinel.
inline bool hasSentinel(const std::int16_t* tokens, std::size_t start, std::size_t k) {
    for (std::size_t j = 0; j < k; ++j) {
        if (tokens[start + j] == kTokenSentinel) return true;
    }
    return false;
}

}  // namespace piru::signal
