// SPDX-License-Identifier: MIT
// GraphStore interface and FlatGraph-backed implementation.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "core/index/flat_graph.hpp"

namespace piru::index {

class GraphStore {
public:
  virtual ~GraphStore() = default;

  virtual std::size_t nodeCount() const = 0;
  virtual std::size_t sequenceLen(std::size_t node_id) const = 0;
  virtual const FlatGraph& flat() const = 0;
};

class FlatGraphStore : public GraphStore {
public:
  FlatGraphStore() = default;
  explicit FlatGraphStore(FlatGraph graph) : graph_(std::move(graph)) {}

  const FlatGraph& flat() const override { return graph_; }

  std::size_t nodeCount() const override { return graph_.nodeCount(); }
  std::size_t sequenceLen(std::size_t node_id) const override {
    return graph_.seqLen(static_cast<std::uint32_t>(node_id));
  }

private:
  FlatGraph graph_;
};

// Legacy alias during migration
using AdjListGraphStore = FlatGraphStore;

}  // namespace piru::index
