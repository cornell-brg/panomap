#pragma once

/**
 * mmap_seed_store.hpp
 *
 * Memory-mapped SeedStore for genome-scale indexes. Instead of loading
 * seed data into heap-allocated vectors, mmaps the pirx file and points
 * directly into it. The OS pages data in/out on demand.
 *
 * Works with any index size -- limited by disk, not RAM. Frequently
 * accessed pages stay in memory (OS page cache), cold pages live on disk.
 *
 * Related:
 *  - seed_store.hpp (FlatSeedStore for small indexes)
 *  - serialization.cpp (pirx format)
 *
 * SPDX-License-Identifier: MIT
 */

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

#include "index/seed_store.hpp"

namespace piru::index {

class MmapSeedStore : public SeedStore {
public:
  // Construct from an already-mmap'd region. The caller owns the mmap lifetime.
  // hashes_ptr: pointer to sorted hash array (n elements)
  // offsets_ptr: pointer to CSR offset array (n+1 elements)
  // entries_ptr: pointer to SeedEntry array (m elements)
  MmapSeedStore(const std::uint64_t* hashes_ptr, std::size_t n_hashes,
                const std::uint32_t* offsets_ptr,
                const SeedEntry* entries_ptr, std::size_t n_entries,
                std::string extractor_name,
                std::map<std::string, std::string> params,
                std::size_t max_hash_frequency,
                std::size_t frequency_threshold,
                double filter_fraction);

  void insert(std::uint64_t, SeedEntry) override {
    // read-only
  }

  SeedHitSpan lookup(std::uint64_t hash) const override;
  std::size_t size() const override { return n_hashes_; }
  std::size_t max_hash_frequency() const override { return max_hash_frequency_; }
  std::size_t frequency_threshold() const override { return frequency_threshold_; }
  void set_frequency_threshold(std::size_t t) override { frequency_threshold_ = t; }
  void recompute_threshold_from_percentile(double percentile) override;
  double filter_fraction() const override { return filter_fraction_; }
  const std::string& extractor_name() const override { return extractor_name_; }
  const std::map<std::string, std::string>& params() const override { return params_; }

private:
  const std::uint64_t* hashes_;
  const std::uint32_t* offsets_;
  const SeedEntry* entries_;
  std::size_t n_hashes_;
  std::size_t n_entries_;

  std::string extractor_name_;
  std::map<std::string, std::string> params_;
  std::size_t max_hash_frequency_;
  std::size_t frequency_threshold_;
  double filter_fraction_;
};

// RAII wrapper for mmap'd file. Holds the mmap and provides byte offsets.
class MmapFile {
public:
  explicit MmapFile(const std::string& path);
  ~MmapFile();

  MmapFile(const MmapFile&) = delete;
  MmapFile& operator=(const MmapFile&) = delete;

  const char* data() const { return data_; }
  std::size_t size() const { return size_; }

  template <typename T>
  const T* at(std::size_t byte_offset) const {
    return reinterpret_cast<const T*>(data_ + byte_offset);
  }

private:
  char* data_{nullptr};
  std::size_t size_{0};
};

}  // namespace piru::index
