/**
 * kmer_hash.hpp
 *
 * Shared k-mer hashing primitives for base-mode indexing and read-side
 * seeding. The indexer (base/index/base_indexer.cpp) and the seeder
 * (base/seeder/base_seeder.cpp) MUST hash identically; centralising the
 * hash function and 2-bit base table here is the single source of truth.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstdint>

namespace piru::core {

// minimap2's hash64 mixer. Spreads small input bits across all output
// bits so lex-min in 2-bit-packed space does not bias minimizer
// selection. `mask` should be (1 << (2*k)) - 1 for k <= 31, or ~0ULL
// for k == 32.
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

// 2-bit encoding: A=0, C=1, G=2, T=3. Anything else (N, lowercase n,
// IUPAC ambiguity codes) returns 4 -- callers treat that as a "reset"
// signal for the rolling k-mer (drop accumulated bases, restart).
// Lowercase a/c/g/t map to the same as uppercase.
inline std::uint8_t base2bit(char b) {
  switch (b) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default:            return 4;
  }
}

}  // namespace piru::core
