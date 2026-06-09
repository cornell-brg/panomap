/**
 * base_seeder.cpp
 *
 * Read-side counterpart to base/index/base_indexer.cpp. Hashing must be
 * bit-identical -- both call core::hash64 with the same kmer_mask, and
 * both encode bases via core::base2bit.
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/seeder/base_seeder.hpp"

#include <algorithm>
#include <limits>

#include "core/util/kmer_hash.hpp"

namespace panomap::base {

namespace {
constexpr std::uint64_t kInvalidHash = std::numeric_limits<std::uint64_t>::max();
}  // namespace

SeedBuffer extract_minimizers(std::string_view bases, const BaseSeederConfig& cfg) {
  SeedBuffer out;
  const std::size_t k = cfg.k;
  const std::size_t w = std::max<std::size_t>(cfg.w, 1);
  if (k == 0 || k > 32 || bases.size() < k) return out;

  const std::uint64_t kmer_mask = (k == 32) ? ~std::uint64_t{0} : ((1ULL << (2 * k)) - 1);
  const std::size_t num_kmers = bases.size() - k + 1;

  /* 1. Roll k-mer over bases; store hash per k-mer start position. */

  std::vector<std::uint64_t> hashes(num_kmers, kInvalidHash);
  std::uint64_t kmer = 0;
  int valid_bases = 0;

  for (std::size_t i = 0; i < bases.size(); ++i) {
    std::uint8_t code = panomap::core::base2bit(bases[i]);
    if (code >= 4) {
      kmer = 0;
      valid_bases = 0;
    } else {
      kmer = ((kmer << 2) | code) & kmer_mask;
      if (valid_bases < static_cast<int>(k)) ++valid_bases;
    }

    if (i + 1 >= k) {
      const std::size_t kmer_idx = i + 1 - k;
      if (valid_bases >= static_cast<int>(k)) {
        hashes[kmer_idx] = panomap::core::hash64(kmer, kmer_mask);
      }
    }
  }

  /* 2. Sliding window of size w, leftmost-min minimizer. Dedup against
   *    last emitted index (consecutive overlapping windows often pick
   *    the same minimum). */

  if (num_kmers < w) return out;

  out.seeds.reserve(num_kmers / w + 1);
  std::size_t last_emitted = std::numeric_limits<std::size_t>::max();
  for (std::size_t i = 0; i + w <= num_kmers; ++i) {
    std::size_t min_idx = i;
    std::uint64_t min_hash = hashes[i];
    for (std::size_t j = i + 1; j < i + w; ++j) {
      if (hashes[j] < min_hash) {
        min_hash = hashes[j];
        min_idx = j;
      }
    }
    if (min_hash == kInvalidHash) continue;
    if (min_idx == last_emitted) continue;
    last_emitted = min_idx;
    out.seeds.push_back(Seed{.hash = min_hash, .position = min_idx, .length = k});
  }

  return out;
}

}  // namespace panomap::base
