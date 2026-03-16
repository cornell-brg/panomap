// SPDX-License-Identifier: MIT
// SeedStore interface and simple hash-map backend.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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

class SeedStore {
public:
  virtual ~SeedStore() = default;

  virtual void insert(std::uint64_t hash, SeedEntry hit) = 0;
  virtual const std::vector<SeedEntry>* lookup(std::uint64_t hash) const = 0;
  virtual std::size_t size() const = 0;
  virtual std::size_t max_hash_frequency() const = 0;
  virtual std::size_t frequency_threshold() const = 0;
  virtual double filter_fraction() const = 0;
  virtual const std::string& extractor_name() const = 0;
  virtual const std::map<std::string, std::string>& params() const = 0;
};

class HashSeedStore : public SeedStore {
public:
  void insert(std::uint64_t hash, SeedEntry hit) override {
    store_[hash].push_back(std::move(hit));
  }

  const std::vector<SeedEntry>* lookup(std::uint64_t hash) const override {
    const auto it = store_.find(hash);
    if (it == store_.end()) return nullptr;
    return &it->second;
  }

  std::size_t size() const override { return store_.size(); }

  std::size_t max_hash_frequency() const override { return max_hash_frequency_; }
  std::size_t frequency_threshold() const override { return frequency_threshold_; }
  double filter_fraction() const override { return filter_fraction_; }
  const std::string& extractor_name() const override { return extractor_name_; }
  const std::map<std::string, std::string>& params() const override { return params_; }

  void set_max_hash_frequency(std::size_t freq) { max_hash_frequency_ = freq; }
  void set_frequency_threshold(std::size_t threshold) { frequency_threshold_ = threshold; }
  void set_filter_fraction(double fraction) { filter_fraction_ = fraction; }

  // Recompute frequency threshold from percentile (0.0-1.0)
  // e.g., 0.5 = 50th percentile, 0.9 = 90th percentile
  void recompute_threshold_from_percentile(double percentile) {
    if (store_.empty()) return;
    std::vector<std::size_t> frequencies;
    frequencies.reserve(store_.size());
    for (const auto& [hash, hits] : store_) {
      frequencies.push_back(hits.size());
    }
    std::sort(frequencies.begin(), frequencies.end());
    double fraction = std::clamp(percentile, 0.0, 1.0);
    std::size_t pos = static_cast<std::size_t>(frequencies.size() * fraction);
    pos = std::min(pos, frequencies.size() - 1);
    frequency_threshold_ = frequencies[pos] + 1;
    filter_fraction_ = percentile;
  }
  void set_extractor_name(std::string name) { extractor_name_ = std::move(name); }
  void set_params(std::map<std::string, std::string> p) { params_ = std::move(p); }

  const std::unordered_map<std::uint64_t, std::vector<SeedEntry>>& data() const { return store_; }
  std::unordered_map<std::uint64_t, std::vector<SeedEntry>>& mutableData() { return store_; }

  // Merge another HashSeedStore into this one.
  // Used for combining thread-local stores after parallel indexing.
  void merge(const HashSeedStore& other) {
    for (const auto& [hash, hits] : other.store_) {
      auto& target = store_[hash];
      target.insert(target.end(), hits.begin(), hits.end());
    }
  }

  // Remove duplicate (node_id, offset) entries from each hit vector.
  // Used after path-guided seeding where shared regions produce duplicates.
  void deduplicate() {
    for (auto& [hash, hits] : store_) {
      if (hits.size() <= 1) continue;
      std::sort(hits.begin(), hits.end());
      hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    }
  }

private:
  std::unordered_map<std::uint64_t, std::vector<SeedEntry>> store_;
  std::size_t max_hash_frequency_{0};
  std::size_t frequency_threshold_{0};
  double filter_fraction_{0.0};
  std::string extractor_name_;
  std::map<std::string, std::string> params_;
};

}  // namespace piru::index
