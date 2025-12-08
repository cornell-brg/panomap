// SPDX-License-Identifier: MIT
// SeedStore interface and simple hash-map backend.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace piru::index {

struct SeedHit {
    std::size_t node_id{0};
    std::size_t offset{0};
};

class SeedStore {
public:
    virtual ~SeedStore() = default;

    virtual void insert(std::uint64_t hash, SeedHit hit) = 0;
    virtual const std::vector<SeedHit>* lookup(std::uint64_t hash) const = 0;
    virtual std::size_t size() const = 0;
    virtual std::size_t max_hash_frequency() const = 0;
    virtual std::size_t frequency_threshold() const = 0;
    virtual uint32_t seed_k() const = 0;
    virtual uint32_t seed_stride() const = 0;
    virtual uint32_t seed_qbits() const = 0;
    virtual double filter_fraction() const = 0;
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
    uint32_t seed_k() const override { return seed_k_; }
    uint32_t seed_stride() const override { return seed_stride_; }
    uint32_t seed_qbits() const override { return seed_qbits_; }
    double filter_fraction() const override { return filter_fraction_; }

    void set_max_hash_frequency(std::size_t freq) { max_hash_frequency_ = freq; }
    void set_frequency_threshold(std::size_t threshold) { frequency_threshold_ = threshold; }
    void set_seed_k(uint32_t k) { seed_k_ = k; }
    void set_seed_stride(uint32_t stride) { seed_stride_ = stride; }
    void set_seed_qbits(uint32_t qbits) { seed_qbits_ = qbits; }
    void set_filter_fraction(double fraction) { filter_fraction_ = fraction; }

    const std::unordered_map<std::uint64_t, std::vector<SeedHit>>& data() const { return store_; }
    std::unordered_map<std::uint64_t, std::vector<SeedHit>>& mutableData() { return store_; }

private:
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> store_;
    std::size_t max_hash_frequency_{0};
    std::size_t frequency_threshold_{0};
    uint32_t seed_k_{0};
    uint32_t seed_stride_{0};
    uint32_t seed_qbits_{0};
    double filter_fraction_{0.0};
};

}  // namespace piru::index
