// SPDX-License-Identifier: MIT
// SeedStore interface and simple hash-map backend.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>

namespace piru::index {

struct SeedHit {
    std::size_t node_id{0};
    std::size_t offset{0};
    std::size_t length{0};     // Seed coverage length (from Seed.length)

    // For sorting and deduplication (ignore length for uniqueness)
    bool operator<(const SeedHit& other) const {
        if (node_id != other.node_id) return node_id < other.node_id;
        return offset < other.offset;
    }
    bool operator==(const SeedHit& other) const {
        return node_id == other.node_id && offset == other.offset;
    }
};

class SeedStore {
public:
    virtual ~SeedStore() = default;

    virtual void insert(std::uint64_t hash, SeedHit hit) = 0;
    virtual const std::vector<SeedHit>* lookup(std::uint64_t hash) const = 0;
    virtual std::size_t size() const = 0;
    virtual std::size_t max_hash_frequency() const = 0;
    virtual std::size_t frequency_threshold() const = 0;
    virtual double filter_fraction() const = 0;
    virtual const std::string& extractor_name() const = 0;
    virtual const std::map<std::string, std::string>& params() const = 0;
};

using SeedStorePtr = std::unique_ptr<SeedStore>;

class HashSeedStore : public SeedStore {
public:
    void insert(std::uint64_t hash, SeedHit hit) override {
        store_[hash].push_back(std::move(hit));
    }

    const std::vector<SeedHit>* lookup(std::uint64_t hash) const override {
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
    void set_extractor_name(std::string name) { extractor_name_ = std::move(name); }
    void set_params(std::map<std::string, std::string> p) { params_ = std::move(p); }

    const std::unordered_map<std::uint64_t, std::vector<SeedHit>>& data() const { return store_; }
    std::unordered_map<std::uint64_t, std::vector<SeedHit>>& mutableData() { return store_; }

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
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> store_;
    std::size_t max_hash_frequency_{0};
    std::size_t frequency_threshold_{0};
    double filter_fraction_{0.0};
    std::string extractor_name_;
    std::map<std::string, std::string> params_;
};

}  // namespace piru::index
