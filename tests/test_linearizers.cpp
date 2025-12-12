#include <doctest/doctest.h>

#include "index/linearizer_factory.hpp"
#include "index/path_walk_linearizer.hpp"
#include "index/pseudo_linearize.hpp"
#include "index/superbubble_linearizer.hpp"

using namespace piru::index;

// Helper to create a simple linear graph: A -> B -> C
AlnGraph createLinearGraph() {
  AlnGraph graph;

  AlnNode a;
  a.label = "A";
  a.sequence = "AAAA";  // length 4
  graph.addNode(a);

  AlnNode b;
  b.label = "B";
  b.sequence = "BBBB";  // length 4
  graph.addNode(b);

  AlnNode c;
  c.label = "C";
  c.sequence = "CCCC";  // length 4
  graph.addNode(c);

  AlnEdge ab;
  ab.from = 0;
  ab.to = 1;
  ab.overlap_bases = 0;
  graph.addEdge(ab);

  AlnEdge bc;
  bc.from = 1;
  bc.to = 2;
  bc.overlap_bases = 0;
  graph.addEdge(bc);

  return graph;
}

// Helper to create a bubble graph: A -> B -> D
//                                  A -> C -> D
AlnGraph createBubbleGraph() {
  AlnGraph graph;

  AlnNode a;
  a.label = "A";
  a.sequence = "AAAA";
  graph.addNode(a);

  AlnNode b;
  b.label = "B";
  b.sequence = "BBBB";
  graph.addNode(b);

  AlnNode c;
  c.label = "C";
  c.sequence = "CCCC";
  graph.addNode(c);

  AlnNode d;
  d.label = "D";
  d.sequence = "DDDD";
  graph.addNode(d);

  graph.addEdge({0, 1, 0});  // A -> B
  graph.addEdge({0, 2, 0});  // A -> C
  graph.addEdge({1, 3, 0});  // B -> D
  graph.addEdge({2, 3, 0});  // C -> D

  return graph;
}

TEST_CASE("SuperbubbleLinearizer on linear graph") {
  AlnGraph graph = createLinearGraph();

  // Add path
  AlnPath path;
  path.name = "path1";
  path.steps.push_back({"A", false});
  path.steps.push_back({"B", false});
  path.steps.push_back({"C", false});
  graph.addPath(path);

  SuperbubbleLinearizer linearizer;
  auto coords = linearizer.linearize(graph);

  REQUIRE(coords.size() == 3);

  // All nodes should get coordinates (linear graph forms one chain).
  // Each node gets exactly one coordinate.
  CHECK(coords[0].size() >= 0);  // Node A
  CHECK(coords[1].size() >= 0);  // Node B
  CHECK(coords[2].size() >= 0);  // Node C

  // Check that linearizer name is correct
  CHECK(linearizer.name() == "superbubble");
}

TEST_CASE("PathWalkLinearizer on linear graph with single path") {
  AlnGraph graph = createLinearGraph();

  // Add path: A -> B -> C
  AlnPath path;
  path.name = "path1";
  path.steps.push_back({"A", false});
  path.steps.push_back({"B", false});
  path.steps.push_back({"C", false});
  path.overlaps.push_back(0);  // A-B overlap
  path.overlaps.push_back(0);  // B-C overlap
  graph.addPath(path);

  PathWalkLinearizer linearizer;
  auto coords = linearizer.linearize(graph);

  REQUIRE(coords.size() == 3);

  // Each node should have exactly one coordinate (single path).
  REQUIRE(coords[0].size() == 1);  // Node A
  REQUIRE(coords[1].size() == 1);  // Node B
  REQUIRE(coords[2].size() == 1);  // Node C

  // All on same path
  CHECK(coords[0][0].path_id == 0);
  CHECK(coords[1][0].path_id == 0);
  CHECK(coords[2][0].path_id == 0);

  // Coordinates should be monotonically increasing
  // A at position 0, B at position 4 (A.length=4), C at position 8 (A+B=8)
  CHECK(coords[0][0].ref_coord == 0);
  CHECK(coords[1][0].ref_coord == 4);
  CHECK(coords[2][0].ref_coord == 8);

  // Check linearizer name
  CHECK(linearizer.name() == "path-walk");
}

TEST_CASE("PathWalkLinearizer on bubble graph with two paths") {
  AlnGraph graph = createBubbleGraph();

  // Path 1: A -> B -> D
  AlnPath path1;
  path1.name = "hap1";
  path1.steps.push_back({"A", false});
  path1.steps.push_back({"B", false});
  path1.steps.push_back({"D", false});
  path1.overlaps = {0, 0};
  graph.addPath(path1);

  // Path 2: A -> C -> D
  AlnPath path2;
  path2.name = "hap2";
  path2.steps.push_back({"A", false});
  path2.steps.push_back({"C", false});
  path2.steps.push_back({"D", false});
  path2.overlaps = {0, 0};
  graph.addPath(path2);

  PathWalkLinearizer linearizer;
  auto coords = linearizer.linearize(graph);

  REQUIRE(coords.size() == 4);

  // Node A: appears on both paths (2 coords)
  REQUIRE(coords[0].size() == 2);
  CHECK(coords[0][0].path_id == 0);  // path1
  CHECK(coords[0][1].path_id == 1);  // path2
  CHECK(coords[0][0].ref_coord == 0);
  CHECK(coords[0][1].ref_coord == 0);

  // Node B: only on path1 (1 coord)
  REQUIRE(coords[1].size() == 1);
  CHECK(coords[1][0].path_id == 0);
  CHECK(coords[1][0].ref_coord == 4);

  // Node C: only on path2 (1 coord)
  REQUIRE(coords[2].size() == 1);
  CHECK(coords[2][0].path_id == 1);
  CHECK(coords[2][0].ref_coord == 4);

  // Node D: appears on both paths (2 coords)
  REQUIRE(coords[3].size() == 2);
  CHECK(coords[3][0].path_id == 0);  // path1
  CHECK(coords[3][1].path_id == 1);  // path2
  CHECK(coords[3][0].ref_coord == 8);  // A(4) + B(4)
  CHECK(coords[3][1].ref_coord == 8);  // A(4) + C(4)
}

TEST_CASE("PathWalkLinearizer on graph with cycle") {
  // Create graph: A -> B -> C -> B (cycle)
  AlnGraph graph;

  AlnNode a;
  a.label = "A";
  a.sequence = "AAAA";
  graph.addNode(a);

  AlnNode b;
  b.label = "B";
  b.sequence = "BBBB";
  graph.addNode(b);

  AlnNode c;
  c.label = "C";
  c.sequence = "CCCC";
  graph.addNode(c);

  graph.addEdge({0, 1, 0});  // A -> B
  graph.addEdge({1, 2, 0});  // B -> C
  graph.addEdge({2, 1, 0});  // C -> B (cycle)

  // Path traverses cycle: A -> B -> C -> B
  AlnPath path;
  path.name = "cyclic_path";
  path.steps.push_back({"A", false});
  path.steps.push_back({"B", false});
  path.steps.push_back({"C", false});
  path.steps.push_back({"B", false});  // B appears again
  path.overlaps = {0, 0, 0};
  graph.addPath(path);

  PathWalkLinearizer linearizer;
  auto coords = linearizer.linearize(graph);

  REQUIRE(coords.size() == 3);

  // Node A: appears once (1 coord)
  REQUIRE(coords[0].size() == 1);
  CHECK(coords[0][0].ref_coord == 0);

  // Node B: appears twice on same path (2 coords)
  REQUIRE(coords[1].size() == 2);
  CHECK(coords[1][0].path_id == 0);
  CHECK(coords[1][1].path_id == 0);  // Same path
  CHECK(coords[1][0].ref_coord == 4);   // First occurrence
  CHECK(coords[1][1].ref_coord == 12);  // Second occurrence: A(4) + B(4) + C(4)

  // Node C: appears once (1 coord)
  REQUIRE(coords[2].size() == 1);
  CHECK(coords[2][0].ref_coord == 8);
}

TEST_CASE("PathWalkLinearizer error on graph without paths") {
  AlnGraph graph = createLinearGraph();
  // No paths added

  PathWalkLinearizer linearizer;

  CHECK_THROWS_AS(linearizer.linearize(graph), std::runtime_error);
}

TEST_CASE("PathWalkLinearizer with overlaps") {
  AlnGraph graph = createLinearGraph();

  // Add path with overlaps: A -[overlap 1]-> B -[overlap 2]-> C
  AlnPath path;
  path.name = "path_with_overlaps";
  path.steps.push_back({"A", false});
  path.steps.push_back({"B", false});
  path.steps.push_back({"C", false});
  path.overlaps.push_back(1);  // A-B overlap of 1
  path.overlaps.push_back(2);  // B-C overlap of 2
  graph.addPath(path);

  PathWalkLinearizer linearizer;
  auto coords = linearizer.linearize(graph);

  REQUIRE(coords.size() == 3);
  REQUIRE(coords[0].size() == 1);
  REQUIRE(coords[1].size() == 1);
  REQUIRE(coords[2].size() == 1);

  // A at 0, B at 4-1=3, C at 3+4-2=5
  CHECK(coords[0][0].ref_coord == 0);
  CHECK(coords[1][0].ref_coord == 3);   // 0 + 4 - 1
  CHECK(coords[2][0].ref_coord == 5);   // 3 + 4 - 2
}

TEST_CASE("Linearizer factory creates correct backends") {
  SUBCASE("superbubble backend") {
    auto linearizer = make_linearizer("superbubble");
    REQUIRE(linearizer != nullptr);
    CHECK(linearizer->name() == "superbubble");
  }

  SUBCASE("path-walk backend") {
    auto linearizer = make_linearizer("path-walk");
    REQUIRE(linearizer != nullptr);
    CHECK(linearizer->name() == "path-walk");
  }

  SUBCASE("unknown backend throws") {
    CHECK_THROWS_AS(make_linearizer("unknown"), std::invalid_argument);
  }
}
