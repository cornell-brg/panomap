// SPDX-License-Identifier: MIT

#include "index/mmap_seed_store.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace piru::index {

// -- MmapFile --

MmapFile::MmapFile(const std::string& path) {
  int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    throw std::runtime_error("Failed to open file for mmap: " + path);
  }

  struct stat st;
  if (fstat(fd, &st) < 0) {
    close(fd);
    throw std::runtime_error("Failed to stat file: " + path);
  }
  size_ = static_cast<std::size_t>(st.st_size);

  data_ = static_cast<char*>(
      mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0));
  close(fd);  // fd no longer needed after mmap

  if (data_ == MAP_FAILED) {
    data_ = nullptr;
    size_ = 0;
    throw std::runtime_error("mmap failed for: " + path);
  }

  // Advise sequential access for initial scan, random for lookups
  madvise(data_, size_, MADV_RANDOM);
}

MmapFile::~MmapFile() {
  if (data_ && size_ > 0) {
    munmap(data_, size_);
  }
}

// -- MmapSeedStore --

MmapSeedStore::MmapSeedStore(const std::uint64_t* hashes_ptr, std::size_t n_hashes,
                               const std::uint32_t* offsets_ptr,
                               const SeedEntry* entries_ptr, std::size_t n_entries,
                               std::string extractor_name,
                               std::map<std::string, std::string> params,
                               std::size_t max_hash_frequency,
                               std::size_t frequency_threshold,
                               double filter_fraction)
    : hashes_(hashes_ptr),
      offsets_(offsets_ptr),
      entries_(entries_ptr),
      n_hashes_(n_hashes),
      n_entries_(n_entries),
      extractor_name_(std::move(extractor_name)),
      params_(std::move(params)),
      max_hash_frequency_(max_hash_frequency),
      frequency_threshold_(frequency_threshold),
      filter_fraction_(filter_fraction) {}

SeedHitSpan MmapSeedStore::lookup(std::uint64_t hash) const {
  // Binary search on sorted hashes (same as FlatSeedStore)
  auto it = std::lower_bound(hashes_, hashes_ + n_hashes_, hash);
  if (it == hashes_ + n_hashes_ || *it != hash) return {};
  auto idx = static_cast<std::size_t>(it - hashes_);
  return {entries_ + offsets_[idx], offsets_[idx + 1] - offsets_[idx]};
}

void MmapSeedStore::recompute_threshold_from_percentile(double percentile) {
  if (n_hashes_ == 0) return;
  std::vector<std::size_t> freqs;
  freqs.reserve(n_hashes_);
  for (std::size_t i = 0; i < n_hashes_; ++i) {
    freqs.push_back(offsets_[i + 1] - offsets_[i]);
  }
  std::sort(freqs.begin(), freqs.end());
  double frac = std::clamp(percentile, 0.0, 1.0);
  auto pos = std::min(static_cast<std::size_t>(freqs.size() * frac), freqs.size() - 1);
  frequency_threshold_ = freqs[pos] + 1;
  filter_fraction_ = percentile;
}

}  // namespace piru::index
