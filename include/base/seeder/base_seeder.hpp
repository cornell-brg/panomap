/**
 * base_seeder.hpp
 *
 * Per-read minimizer extraction for base mode. Mirrors the indexer:
 * directional 2-bit k-mer rolling, minimap2 hash64, sliding window of w
 * with leftmost tie-break. Hashes the read in basecalled order only --
 * the directional FlatGraph already covers both strand orientations, so
 * there's no read-side reverse complement (matches the signal-mode
 * contract: index handles RC, read does not).
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace panomap::base {

// Same shape as signal::Seed; kept as a base-namespace type to honour the
// "no signal:: dependency in panomap_base" architectural rule.
struct Seed {
  std::uint64_t hash{0};
  std::size_t position{0};
  std::size_t length{0};
};

struct SeedBuffer {
  std::vector<Seed> seeds;
};

struct BaseSeederConfig {
  std::size_t k{15};
  std::size_t w{10};
};

// Extract minimizers from `bases` (ASCII, may contain N or lowercase).
// On N: rolls reset, the affected k-mer windows produce no minimizers.
// Output positions are k-mer start offsets in the read; length = k.
SeedBuffer extract_minimizers(std::string_view bases, const BaseSeederConfig& cfg);

}  // namespace panomap::base
