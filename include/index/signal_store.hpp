// SPDX-License-Identifier: MIT
// SignalStore interface and simple in-memory backend.

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "signal/signal_types.hpp"

namespace piru::index {

class SignalStore {
public:
    virtual ~SignalStore() = default;

    virtual std::size_t size() const = 0;
    virtual const piru::signal::AlignmentQuantizedSignal* get(std::size_t node_id) const = 0;
};

using SignalStorePtr = std::unique_ptr<SignalStore>;

class VectorSignalStore : public SignalStore {
public:
    VectorSignalStore() = default;
    explicit VectorSignalStore(std::vector<piru::signal::AlignmentQuantizedSignal> signals)
        : signals_(std::move(signals)) {}

    std::size_t size() const override { return signals_.size(); }

    const piru::signal::AlignmentQuantizedSignal* get(std::size_t node_id) const override {
        if (node_id >= signals_.size()) return nullptr;
        return &signals_[node_id];
    }

    std::vector<piru::signal::AlignmentQuantizedSignal>& mutableSignals() { return signals_; }

private:
    std::vector<piru::signal::AlignmentQuantizedSignal> signals_;
};

}  // namespace piru::index
