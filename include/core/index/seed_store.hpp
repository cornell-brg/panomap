// SPDX-License-Identifier: MIT
// SeedStore interface. Canonical implementation: BucketSeedStore (bucket_seed_store.hpp).

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

namespace piru::index {

struct SeedEntry {
  std::uint32_t node_id{0};
  std::uint32_t offset{0};

  bool operator<(const SeedEntry& other) const {
    if (node_id != other.node_id) return node_id < other.node_id;
    return offset < other.offset;
  }
  bool operator==(const SeedEntry& other) const {
    return node_id == other.node_id && offset == other.offset;
  }
};

// Lightweight view into a contiguous array of SeedEntry.
struct SeedHitSpan {
  const SeedEntry* data{nullptr};
  std::size_t count{0};
};

class SeedStore {
public:
  virtual ~SeedStore() = default;

  virtual void insert(std::uint64_t hash, SeedEntry hit) = 0;
  virtual SeedHitSpan lookup(std::uint64_t hash) const = 0;
  virtual std::size_t size() const = 0;
  virtual std::size_t max_hash_frequency() const = 0;
  virtual std::size_t frequency_threshold() const = 0;
  virtual void set_frequency_threshold(std::size_t threshold) = 0;
  virtual void recompute_threshold_from_percentile(double percentile) = 0;

  // Mirror minimap2's mm_idx_cal_max_occ: threshold = (1 - top_frac)
  // percentile of frequencies + 1, then clamped to [min_occ, max_occ].
  // Default minimap2 args: top_frac=2e-4, min=10, max=1e6 -> ~694 on hg38
  // (vs ~45 from p0.99 percentile, which over-filters at human scale).
  virtual void recompute_threshold_from_top_frac(double top_frac, std::size_t min_occ,
                                                 std::size_t max_occ) = 0;

  virtual double filter_fraction() const = 0;
  virtual const std::string& extractor_name() const = 0;
  virtual const std::map<std::string, std::string>& params() const = 0;
};

}  // namespace piru::index
