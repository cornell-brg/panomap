// SPDX-License-Identifier: MIT

#include "alignment/chain_aligner.hpp"
#include "alignment/segment_aligner.hpp"
#include "alignment/signal_utils.hpp"

#include <doctest/doctest.h>

using namespace piru::alignment;
using namespace piru::signal;

// =============================================================================
// Signal Utilities Tests
// =============================================================================

TEST_CASE("signalLength returns correct size for float32") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kFloat32;
  sig.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f};

  CHECK(signalLength(sig) == 5);
}

TEST_CASE("signalLength returns correct size for int16") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kInt16;
  sig.data = std::vector<int16_t>{100, 200, 300};

  CHECK(signalLength(sig) == 3);
}

TEST_CASE("sliceSignal extracts correct range for float32") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kFloat32;
  sig.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f};

  auto sliced = sliceSignal(sig, 1, 4);

  CHECK(sliced.kind == AlignmentQuantizationKind::kFloat32);
  auto* vec = std::get_if<std::vector<float>>(&sliced.data);
  REQUIRE(vec != nullptr);
  REQUIRE(vec->size() == 3);
  CHECK((*vec)[0] == 2.0f);
  CHECK((*vec)[1] == 3.0f);
  CHECK((*vec)[2] == 4.0f);
}

TEST_CASE("sliceSignal handles empty range") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kFloat32;
  sig.data = std::vector<float>{1.0f, 2.0f, 3.0f};

  auto sliced = sliceSignal(sig, 2, 2);

  CHECK(signalLength(sliced) == 0);
}

TEST_CASE("sliceSignal handles out-of-bounds gracefully") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kFloat32;
  sig.data = std::vector<float>{1.0f, 2.0f, 3.0f};

  auto sliced = sliceSignal(sig, 1, 100);

  auto* vec = std::get_if<std::vector<float>>(&sliced.data);
  REQUIRE(vec != nullptr);
  CHECK(vec->size() == 2);  // Only indices 1, 2
}

TEST_CASE("concatSignals combines two signals") {
  AlignmentQuantizedSignal a, b;
  a.kind = AlignmentQuantizationKind::kFloat32;
  a.data = std::vector<float>{1.0f, 2.0f};
  b.kind = AlignmentQuantizationKind::kFloat32;
  b.data = std::vector<float>{3.0f, 4.0f, 5.0f};

  auto combined = concatSignals(a, b);

  CHECK(combined.kind == AlignmentQuantizationKind::kFloat32);
  auto* vec = std::get_if<std::vector<float>>(&combined.data);
  REQUIRE(vec != nullptr);
  REQUIRE(vec->size() == 5);
  CHECK((*vec)[0] == 1.0f);
  CHECK((*vec)[1] == 2.0f);
  CHECK((*vec)[2] == 3.0f);
  CHECK((*vec)[3] == 4.0f);
  CHECK((*vec)[4] == 5.0f);
}

TEST_CASE("concatSignals throws on mismatched kinds") {
  AlignmentQuantizedSignal a, b;
  a.kind = AlignmentQuantizationKind::kFloat32;
  a.data = std::vector<float>{1.0f};
  b.kind = AlignmentQuantizationKind::kInt16;
  b.data = std::vector<int16_t>{100};

  CHECK_THROWS_AS(concatSignals(a, b), std::invalid_argument);
}

TEST_CASE("signalValueAt returns correct float value") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kInt16;
  sig.data = std::vector<int16_t>{100, -200, 300};

  CHECK(signalValueAt(sig, 0) == 100.0f);
  CHECK(signalValueAt(sig, 1) == -200.0f);
  CHECK(signalValueAt(sig, 2) == 300.0f);
}

TEST_CASE("signalValueAt throws on out-of-bounds") {
  AlignmentQuantizedSignal sig;
  sig.kind = AlignmentQuantizationKind::kFloat32;
  sig.data = std::vector<float>{1.0f, 2.0f};

  CHECK_THROWS_AS(signalValueAt(sig, 5), std::out_of_range);
}

TEST_CASE("sameQuantizationKind detects matching kinds") {
  AlignmentQuantizedSignal a, b, c;
  a.kind = AlignmentQuantizationKind::kFloat32;
  a.data = std::vector<float>{};
  b.kind = AlignmentQuantizationKind::kFloat32;
  b.data = std::vector<float>{};
  c.kind = AlignmentQuantizationKind::kInt16;
  c.data = std::vector<int16_t>{};

  CHECK(sameQuantizationKind(a, b));
  CHECK_FALSE(sameQuantizationKind(a, c));
}

// =============================================================================
// Core Types Tests
// =============================================================================

TEST_CASE("GraphPosition equality") {
  GraphPosition a{10, 5};
  GraphPosition b{10, 5};
  GraphPosition c{10, 6};
  GraphPosition d{11, 5};

  CHECK(a == b);
  CHECK_FALSE(a == c);
  CHECK_FALSE(a == d);
}

TEST_CASE("DtwResult valid check") {
  DtwResult valid_result;
  valid_result.cost = 10.0f;
  CHECK(valid_result.valid());

  DtwResult invalid_result;
  invalid_result.cost = std::numeric_limits<float>::infinity();
  CHECK_FALSE(invalid_result.valid());
}

// =============================================================================
// PathGuidedDtwAligner Tests
// =============================================================================

namespace {

// Simple mock graph for testing
class MockGraphStore : public piru::index::GraphStore {
 public:
  struct Node {
    std::string sequence;
    std::vector<size_t> outgoing;
    std::vector<size_t> incoming;
    std::optional<int64_t> chain_id;
    std::optional<int64_t> linear_pos;
  };

  std::vector<Node> nodes;

  size_t nodeCount() const override { return nodes.size(); }
  const std::string& sequence(size_t id) const override { return nodes[id].sequence; }
  const std::vector<size_t>& outgoing(size_t id) const override { return nodes[id].outgoing; }
  const std::vector<size_t>& incoming(size_t id) const override { return nodes[id].incoming; }
  std::optional<int64_t> chainId(size_t id) const override { return nodes[id].chain_id; }
  std::optional<int64_t> linearPosition(size_t id) const override { return nodes[id].linear_pos; }
};

// Simple mock signal store for testing
class MockSignalStore : public piru::index::SignalStore {
 public:
  std::vector<AlignmentQuantizedSignal> signals;

  size_t size() const override { return signals.size(); }
  const AlignmentQuantizedSignal* get(size_t id) const override {
    if (id >= signals.size()) return nullptr;
    return &signals[id];
  }
};

}  // namespace

TEST_CASE("PathGuidedDtwAligner with simple linear graph") {
  // Create a simple linear graph: node0 -> node1 -> node2
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AAA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BBB", .outgoing = {2}, .incoming = {0}, .chain_id = 0, .linear_pos = 3},
      {.sequence = "CCC", .outgoing = {}, .incoming = {1}, .chain_id = 0, .linear_pos = 6},
  };

  // Create signals for each node (float32)
  MockSignalStore signals;
  signals.signals.resize(3);

  // Node 0: [1.0, 2.0, 3.0]
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f, 3.0f};

  // Node 1: [4.0, 5.0, 6.0]
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{4.0f, 5.0f, 6.0f};

  // Node 2: [7.0, 8.0, 9.0]
  signals.signals[2].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[2].data = std::vector<float>{7.0f, 8.0f, 9.0f};

  // Create query signal that matches part of the path
  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{2.0f, 3.0f, 4.0f, 5.0f};  // Should align to node0[1:] + node1[:2]

  // Create aligner
  auto aligner = makePathGuidedDtwAligner();
  REQUIRE(aligner != nullptr);
  CHECK(aligner->name() == "PathGuidedDtwAligner");

  // Align from node0 offset 1 to node1 offset 2
  Anchor start{.graph_pos = {0, 1}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 4};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost < 1.0f);  // Should be near-perfect match
  CHECK(result.reached_target);
  CHECK_FALSE(result.path.empty());
}

TEST_CASE("PathGuidedDtwAligner with identical signals") {
  // Simple graph: node0 -> node1
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 2},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};

  // Query that exactly matches ref
  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  auto aligner = makePathGuidedDtwAligner();
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 4};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));  // Exact match
  CHECK(result.reached_target);
  CHECK(result.path.size() == 4);  // One position per signal sample
}

TEST_CASE("PathGuidedDtwAligner with int16 signals") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "A", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "B", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 1},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kInt16;
  signals.signals[0].data = std::vector<int16_t>{100, 200};
  signals.signals[1].kind = AlignmentQuantizationKind::kInt16;
  signals.signals[1].data = std::vector<int16_t>{300, 400};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kInt16;
  query.data = std::vector<int16_t>{100, 200, 300, 400};

  auto aligner = makePathGuidedDtwAligner();
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 4};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));
}

TEST_CASE("PathGuidedDtwAligner with single node") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AAAA", .outgoing = {}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
  };

  MockSignalStore signals;
  signals.signals.resize(1);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{2.0f, 3.0f};

  auto aligner = makePathGuidedDtwAligner();
  Anchor start{.graph_pos = {0, 1}, .query_pos = 0};
  Anchor target{.graph_pos = {0, 3}, .query_pos = 2};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));
}

TEST_CASE("PathGuidedDtwAligner with no path returns invalid") {
  MockGraphStore graph;
  // Two disconnected nodes
  graph.nodes = {
      {.sequence = "AA", .outgoing = {}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {}, .incoming = {}, .chain_id = 1, .linear_pos = 0},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f};

  auto aligner = makePathGuidedDtwAligner();
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 2};  // Node 1 not reachable from node 0

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK_FALSE(result.valid());
  CHECK_FALSE(result.reached_target);
}

// =============================================================================
// RadiusGdtwAligner Tests
// =============================================================================

TEST_CASE("RadiusGdtwAligner with simple linear graph") {
  // Create a simple linear graph: node0 -> node1 -> node2
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AAA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BBB", .outgoing = {2}, .incoming = {0}, .chain_id = 0, .linear_pos = 3},
      {.sequence = "CCC", .outgoing = {}, .incoming = {1}, .chain_id = 0, .linear_pos = 6},
  };

  MockSignalStore signals;
  signals.signals.resize(3);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f, 3.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{4.0f, 5.0f, 6.0f};
  signals.signals[2].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[2].data = std::vector<float>{7.0f, 8.0f, 9.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{2.0f, 3.0f, 4.0f, 5.0f};

  auto aligner = makeRadiusGdtwAligner(10);
  REQUIRE(aligner != nullptr);
  CHECK(aligner->name() == "RadiusGdtwAligner");

  Anchor start{.graph_pos = {0, 1}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 4};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost < 1.0f);  // Should be near-perfect match
  CHECK_FALSE(result.path.empty());
}

TEST_CASE("RadiusGdtwAligner with identical signals") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 2},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  auto aligner = makeRadiusGdtwAligner(10);
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 2}, .query_pos = 4};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));  // Exact match
  CHECK(result.path.size() == 4);  // One position per signal sample
}

TEST_CASE("RadiusGdtwAligner with branching graph") {
  // Graph with branch: node0 -> node1, node0 -> node2
  // Query should find best path
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "A", .outgoing = {1, 2}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "B", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 1},
      {.sequence = "C", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 1},
  };

  MockSignalStore signals;
  signals.signals.resize(3);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{2.0f};  // Branch 1
  signals.signals[2].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[2].data = std::vector<float>{10.0f};  // Branch 2 (different value)

  // Query matches branch 1 better
  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f};

  auto aligner = makeRadiusGdtwAligner(10);
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {1, 1}, .query_pos = 2};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));  // Exact match on branch 1
  REQUIRE(result.path.size() == 2);
  CHECK(result.path[0].node_id == 0);
  CHECK(result.path[1].node_id == 1);  // Should take branch 1
}

TEST_CASE("RadiusGdtwAligner with single node") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AAAA", .outgoing = {}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
  };

  MockSignalStore signals;
  signals.signals.resize(1);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{2.0f, 3.0f};

  auto aligner = makeRadiusGdtwAligner(10);
  Anchor start{.graph_pos = {0, 1}, .query_pos = 0};
  Anchor target{.graph_pos = {0, 3}, .query_pos = 2};

  auto result = aligner->align(graph, signals, query, start, target);

  CHECK(result.valid());
  CHECK(result.cost == doctest::Approx(0.0f));
}

TEST_CASE("RadiusGdtwAligner respects radius limit") {
  // Long chain: node0 -> node1 -> node2 -> node3 -> node4
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "A", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "B", .outgoing = {2}, .incoming = {0}, .chain_id = 0, .linear_pos = 1},
      {.sequence = "C", .outgoing = {3}, .incoming = {1}, .chain_id = 0, .linear_pos = 2},
      {.sequence = "D", .outgoing = {4}, .incoming = {2}, .chain_id = 0, .linear_pos = 3},
      {.sequence = "E", .outgoing = {}, .incoming = {3}, .chain_id = 0, .linear_pos = 4},
  };

  MockSignalStore signals;
  signals.signals.resize(5);
  for (size_t i = 0; i < 5; ++i) {
    signals.signals[i].kind = AlignmentQuantizationKind::kFloat32;
    signals.signals[i].data = std::vector<float>{static_cast<float>(i + 1)};
  }

  // Short query
  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f};

  // Use small radius_buffer=0, so radius = query_len + 0 = 2
  auto aligner = makeRadiusGdtwAligner(0);
  Anchor start{.graph_pos = {0, 0}, .query_pos = 0};
  Anchor target{.graph_pos = {4, 1}, .query_pos = 2};  // Target far away

  auto result = aligner->align(graph, signals, query, start, target);

  // Should still produce valid alignment within radius
  CHECK(result.valid());
  // Target (node4) is beyond radius, so reached_target should be false
  CHECK_FALSE(result.reached_target);
}

// =============================================================================
// ChainAligner Tests
// =============================================================================

TEST_CASE("ChainAligner with two anchors (single segment)") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 2},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0},
      {.graph_pos = {1, 2}, .query_pos = 4}
  };

  ChainAligner aligner;
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK(result.valid());
  CHECK(result.total_cost == doctest::Approx(0.0f));
  CHECK(result.segments_aligned == 1);
  CHECK(result.segments_reached_target == 1);
  CHECK(result.path.size() == 4);
}

TEST_CASE("ChainAligner with three anchors (two segments)") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {2}, .incoming = {0}, .chain_id = 0, .linear_pos = 2},
      {.sequence = "CC", .outgoing = {}, .incoming = {1}, .chain_id = 0, .linear_pos = 4},
  };

  MockSignalStore signals;
  signals.signals.resize(3);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};
  signals.signals[2].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[2].data = std::vector<float>{5.0f, 6.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};

  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0},
      {.graph_pos = {1, 0}, .query_pos = 2},
      {.graph_pos = {2, 2}, .query_pos = 6}
  };

  ChainAligner aligner;
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK(result.valid());
  CHECK(result.total_cost == doctest::Approx(0.0f));
  CHECK(result.segments_aligned == 2);
}

TEST_CASE("ChainAligner with single anchor returns invalid") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
  };

  MockSignalStore signals;
  signals.signals.resize(1);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f};

  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0}
  };

  ChainAligner aligner;
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK_FALSE(result.valid());
}

TEST_CASE("ChainAligner with radius backend") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BB", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 2},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{3.0f, 4.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};

  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0},
      {.graph_pos = {1, 2}, .query_pos = 4}
  };

  ChainAlignerConfig config;
  config.backend = AlignerBackend::kRadius;
  config.radius_buffer = 10;

  ChainAligner aligner(config);
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK(result.valid());
  CHECK(result.total_cost == doctest::Approx(0.0f));
}

TEST_CASE("ChainAligner auto backend selects appropriately") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AAAA", .outgoing = {1}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
      {.sequence = "BBBB", .outgoing = {}, .incoming = {0}, .chain_id = 0, .linear_pos = 4},
  };

  MockSignalStore signals;
  signals.signals.resize(2);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f};
  signals.signals[1].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[1].data = std::vector<float>{5.0f, 6.0f, 7.0f, 8.0f};

  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};

  // Three anchors: first two on same node (PathGuided), last crosses nodes (Radius)
  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0},
      {.graph_pos = {0, 2}, .query_pos = 2},  // Same node as first
      {.graph_pos = {1, 4}, .query_pos = 8}   // Different node
  };

  ChainAlignerConfig config;
  config.backend = AlignerBackend::kAuto;
  config.radius_buffer = 10;

  ChainAligner aligner(config);
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK(result.valid());
  CHECK(result.segments_aligned == 2);
}

TEST_CASE("ChainAligner normalizedCost calculation") {
  MockGraphStore graph;
  graph.nodes = {
      {.sequence = "AA", .outgoing = {}, .incoming = {}, .chain_id = 0, .linear_pos = 0},
  };

  MockSignalStore signals;
  signals.signals.resize(1);
  signals.signals[0].kind = AlignmentQuantizationKind::kFloat32;
  signals.signals[0].data = std::vector<float>{1.0f, 2.0f};

  // Query slightly off from reference
  AlignmentQuantizedSignal query;
  query.kind = AlignmentQuantizationKind::kFloat32;
  query.data = std::vector<float>{1.5f, 2.5f};

  std::vector<Anchor> anchors = {
      {.graph_pos = {0, 0}, .query_pos = 0},
      {.graph_pos = {0, 2}, .query_pos = 2}
  };

  ChainAligner aligner;
  auto result = aligner.align(graph, signals, query, anchors);

  CHECK(result.valid());
  CHECK(result.total_cost > 0.0f);

  float normalized = result.normalizedCost(2);
  CHECK(normalized == doctest::Approx(result.total_cost / 2.0f));
}
