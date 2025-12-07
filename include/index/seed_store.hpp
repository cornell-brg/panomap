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

    const std::unordered_map<std::uint64_t, std::vector<SeedHit>>& data() const { return store_; }
    std::unordered_map<std::uint64_t, std::vector<SeedHit>>& mutableData() { return store_; }

private:
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> store_;
};

}  // namespace piru::index
