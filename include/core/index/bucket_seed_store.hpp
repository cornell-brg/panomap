#pragma once

/**
 * bucket_seed_store.hpp
 *
 * Bucket-finalized seed store. Buckets are the final index structure --
 * no global merge or CSR at runtime.
 *
 * Each bucket contains:
 *   - keys[]:    sorted unique hash values
 *   - counts[]:  hit count per hash
 *   - offsets[]: offset into entries[] per hash
 *   - entries[]: flat array of all SeedEntry hits (stable memory)
 *
 * Lookup: pick bucket from hash low bits, binary search within bucket.
 *
 * .pirx serialization writes bucket-native format directly.
 *
 * Related:
 *  - seed_store.hpp (SeedStore base interface)
 *  - bucket_indexer.cpp (builds this store)
 *
 * SPDX-License-Identifier: MIT
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "core/index/seed_store.hpp"

namespace piru::index {

/* Per-bucket finalized data. After build, each bucket is read-only.
 *
 * Layout (matches RH2 ri_idx_bucket_s concept):
 *   keys[i]        -- unique hash values in this bucket, sorted
 *   counts[i]      -- hit count for this hash (1 = singleton)
 *   offsets[i]     -- offset into entries[] for this hash's hits
 *   entries[]      -- flat array of all SeedEntry (singletons + multi-hit)
 *
 * All hits live in stable memory (entries[]). No bit-packing, no
 * thread_local. Singletons are just entries[offset] with count=1.
 */
struct Bucket {
  std::vector<std::uint64_t> keys;     // sorted unique hashes
  std::vector<std::uint32_t> counts;   // hit count per hash
  std::vector<std::uint32_t> offsets;  // offset into entries per hash
  std::vector<SeedEntry> entries;      // all hits, stable memory

  /* Lookup a hash in this bucket. Returns span of matching entries. */
  SeedHitSpan lookup(std::uint64_t hash) const {
    auto it = std::lower_bound(keys.begin(), keys.end(), hash);
    if (it == keys.end() || *it != hash) return {};
    auto idx = static_cast<std::size_t>(it - keys.begin());
    return {entries.data() + offsets[idx], counts[idx]};
  }
};

class BucketSeedStore : public SeedStore {
public:
  BucketSeedStore(std::vector<Bucket> buckets, std::uint32_t bucket_bits,
                  std::string extractor_name, std::map<std::string, std::string> params,
                  std::size_t max_hash_frequency, std::size_t frequency_threshold,
                  double filter_fraction)
      : buckets_(std::move(buckets)),
        bucket_bits_(bucket_bits),
        bucket_mask_((1ULL << bucket_bits) - 1),
        extractor_name_(std::move(extractor_name)),
        params_(std::move(params)),
        max_hash_frequency_(max_hash_frequency),
        frequency_threshold_(frequency_threshold),
        filter_fraction_(filter_fraction) {}

  void insert(std::uint64_t, SeedEntry) override {
    throw std::logic_error("BucketSeedStore is read-only");
  }

  SeedHitSpan lookup(std::uint64_t hash) const override {
    std::size_t bi = hash & bucket_mask_;
    return buckets_[bi].lookup(hash);
  }

  std::size_t size() const override {
    std::size_t total = 0;
    for (const auto& b : buckets_) total += b.keys.size();
    return total;
  }

  std::size_t max_hash_frequency() const override { return max_hash_frequency_; }
  std::size_t frequency_threshold() const override { return frequency_threshold_; }
  void set_frequency_threshold(std::size_t threshold) override { frequency_threshold_ = threshold; }

  void recompute_threshold_from_percentile(double percentile) override {
    std::vector<std::size_t> freqs;
    for (const auto& b : buckets_) {
      for (std::size_t i = 0; i < b.keys.size(); ++i) {
        freqs.push_back(b.counts[i]);
      }
    }
    if (freqs.empty()) return;
    std::sort(freqs.begin(), freqs.end());
    double frac = std::clamp(percentile, 0.0, 1.0);
    auto pos = std::min(static_cast<std::size_t>(freqs.size() * frac), freqs.size() - 1);
    frequency_threshold_ = freqs[pos] + 1;
    filter_fraction_ = percentile;
  }

  // Minimap2-style mid_occ. top_frac in [0,1]. Default 2e-4 keeps all but
  // the top 0.02% most frequent seeds. Threshold then clamped to
  // [min_occ, max_occ]. See related/minimap2/index.c:mm_idx_cal_max_occ.
  void recompute_threshold_from_top_frac(double top_frac, std::size_t min_occ,
                                         std::size_t max_occ) override {
    if (top_frac <= 0.0) {
      frequency_threshold_ = std::numeric_limits<std::size_t>::max();
      filter_fraction_ = 0.0;
      return;
    }
    std::vector<std::size_t> freqs;
    for (const auto& b : buckets_) {
      for (std::size_t i = 0; i < b.keys.size(); ++i) {
        freqs.push_back(b.counts[i]);
      }
    }
    if (freqs.empty()) return;
    std::sort(freqs.begin(), freqs.end());
    // (1 - top_frac) percentile = the boundary above which the top top_frac
    // most frequent seeds live. mm_idx_cal_max_occ adds +1 so the threshold
    // is exclusive.
    double clamped = std::clamp(1.0 - top_frac, 0.0, 1.0);
    auto pos = std::min(static_cast<std::size_t>(freqs.size() * clamped), freqs.size() - 1);
    std::size_t thres = freqs[pos] + 1;
    if (thres < min_occ) thres = min_occ;
    if (max_occ > min_occ && thres > max_occ) thres = max_occ;
    frequency_threshold_ = thres;
    filter_fraction_ = top_frac;
  }

  double filter_fraction() const override { return filter_fraction_; }
  const std::string& extractor_name() const override { return extractor_name_; }
  const std::map<std::string, std::string>& params() const override { return params_; }

  // Raw access for serialization
  std::uint32_t bucket_bits() const { return bucket_bits_; }
  std::size_t num_buckets() const { return buckets_.size(); }
  const Bucket& bucket(std::size_t i) const { return buckets_[i]; }

  /* Build a finalized bucket from a sorted array of (hash, entry) pairs.
   * The input must be sorted by hash. Called once per bucket during indexing. */
  static Bucket finalize_bucket(
      const std::vector<std::pair<std::uint64_t, SeedEntry>>& sorted_pairs) {
    Bucket b;
    if (sorted_pairs.empty()) return b;

    // Count unique hashes
    std::size_t n_unique = 0;
    for (std::size_t i = 0; i < sorted_pairs.size();) {
      std::size_t j = i + 1;
      while (j < sorted_pairs.size() && sorted_pairs[j].first == sorted_pairs[i].first) ++j;
      ++n_unique;
      i = j;
    }

    b.keys.reserve(n_unique);
    b.counts.reserve(n_unique);
    b.offsets.reserve(n_unique);
    b.entries.reserve(sorted_pairs.size());

    for (std::size_t i = 0; i < sorted_pairs.size();) {
      std::uint64_t hash = sorted_pairs[i].first;
      std::size_t j = i + 1;
      while (j < sorted_pairs.size() && sorted_pairs[j].first == hash) ++j;

      b.keys.push_back(hash);
      b.counts.push_back(static_cast<std::uint32_t>(j - i));
      b.offsets.push_back(static_cast<std::uint32_t>(b.entries.size()));
      for (std::size_t k = i; k < j; ++k) {
        b.entries.push_back(sorted_pairs[k].second);
      }
      i = j;
    }

    return b;
  }

private:
  std::vector<Bucket> buckets_;
  std::uint32_t bucket_bits_;
  std::uint64_t bucket_mask_;
  std::string extractor_name_;
  std::map<std::string, std::string> params_;
  std::size_t max_hash_frequency_{0};
  std::size_t frequency_threshold_{0};
  double filter_fraction_{0.0};
};

}  // namespace piru::index
