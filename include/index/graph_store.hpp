// SPDX-License-Identifier: MIT
// GraphStore interface and simple adjacency-list backend.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"

namespace piru::index {

class GraphStore {
public:
  virtual ~GraphStore() = default;

  virtual std::size_t nodeCount() const = 0;
  virtual const std::string& sequence(std::size_t node_id) const = 0;
  virtual const std::vector<std::size_t>& outgoing(std::size_t node_id) const = 0;
  virtual const std::vector<std::size_t>& incoming(std::size_t node_id) const = 0;
};

class AdjListGraphStore : public GraphStore {
public:
  AdjListGraphStore() = default;
  explicit AdjListGraphStore(AlnGraph graph) : graph_(std::move(graph)) {}

  const AlnGraph& graph() const { return graph_; }
  AlnGraph& mutableGraph() { return graph_; }

  std::size_t nodeCount() const override { return graph_.nodeCount(); }
  const std::string& sequence(std::size_t node_id) const override {
    return graph_.node(node_id).sequence;
  }
  const std::vector<std::size_t>& outgoing(std::size_t node_id) const override {
    return graph_.outgoing(node_id);
  }
  const std::vector<std::size_t>& incoming(std::size_t node_id) const override {
    return graph_.incoming(node_id);
  }

private:
  AlnGraph graph_;
};

}  // namespace piru::index
