#include <doctest/doctest.h>

#include "index/path_guided_transform.hpp"

using namespace piru::index;
using namespace piru::io;

//------------------------------------------------------------------------------
// Test 1: Single linear path, verify k-1 context added correctly
//------------------------------------------------------------------------------
TEST_CASE("PathGuidedTransform: single linear path with k-1 context") {
  ImportedGraph imported;
  imported.flavor = ImportedGraphFlavor::kVg;

  // Create nodes: 1 -> 2 -> 3
  imported.add_node({"1", "AAAA"});
  imported.add_node({"2", "CCCC"});
  imported.add_node({"3", "GGGG"});

  // Add edges
  imported.add_edge({"1", "2", false, false, "0", std::nullopt});
  imported.add_edge({"2", "3", false, false, "0", std::nullopt});

  // Add path covering all nodes
  ImportedPath path;
  path.name = "path1";
  path.steps.push_back({"1", false});
  path.steps.push_back({"2", false});
  path.steps.push_back({"3", false});
  imported.add_path(path);

  // Transform with pore_k = 5 (k-1 = 4)
  PathGuidedTransform transform;
  AlnGraph result = transform.apply(imported, 0, 5);

  // Verify paths were created
  CHECK(result.pathCount() == 1);

  // Verify coverage statistics
  auto stats = transform.getStats();
  CHECK(stats.original_node_count == 3);
  CHECK(stats.uncovered_node_count == 0);  // All nodes covered by path

  // Verify nodes have k-1 context added
  // Node 1 should have context from node 2 (first 4 bases of "CCCC")
  // Node 2 should have context from node 3 (first 4 bases of "GGGG")
  // Node 3 has no successor, so no context added
  bool found_node1_with_context = false;
  bool found_node2_with_context = false;

  for (std::size_t i = 0; i < result.nodeCount(); ++i) {
    const AlnNode& node = result.node(i);
    if (node.original_id == "1" && !node.is_reverse) {
      if (node.sequence == "AAAACCCC") {  // AAAA + CCCC
        found_node1_with_context = true;
      }
    }
    if (node.original_id == "2" && !node.is_reverse) {
      if (node.sequence == "CCCCGGGG") {  // CCCC + GGGG
        found_node2_with_context = true;
      }
    }
  }

  CHECK(found_node1_with_context);
  CHECK(found_node2_with_context);
}

//------------------------------------------------------------------------------
// Test 2: Two paths sharing node with same successor → node reused
//------------------------------------------------------------------------------
TEST_CASE("PathGuidedTransform: two paths sharing node with same context") {
  ImportedGraph imported;
  imported.flavor = ImportedGraphFlavor::kVg;

  // Linear graph: 1 -> 2 -> 3
  imported.add_node({"1", "AAAA"});
  imported.add_node({"2", "CCCC"});
  imported.add_node({"3", "GGGG"});

  imported.add_edge({"1", "2", false, false, "0", std::nullopt});
  imported.add_edge({"2", "3", false, false, "0", std::nullopt});

  // Path 1: 1 -> 2 -> 3
  ImportedPath path1;
  path1.name = "path1";
  path1.steps.push_back({"1", false});
  path1.steps.push_back({"2", false});
  path1.steps.push_back({"3", false});
  imported.add_path(path1);

  // Path 2: 1 -> 2 (same nodes, same successors)
  ImportedPath path2;
  path2.name = "path2";
  path2.steps.push_back({"1", false});
  path2.steps.push_back({"2", false});
  imported.add_path(path2);

  PathGuidedTransform transform;
  AlnGraph result = transform.apply(imported, 0, 5);

  // Both paths should share the same node variants (same context)
  CHECK(result.pathCount() == 2);

  auto stats = transform.getStats();
  CHECK(stats.original_node_count == 3);
  CHECK(stats.uncovered_node_count == 0);
}

//------------------------------------------------------------------------------
// Test 3: Two paths sharing node with different successors → node duplicated
//------------------------------------------------------------------------------
TEST_CASE("PathGuidedTransform: two paths with different successor contexts") {
  ImportedGraph imported;
  imported.flavor = ImportedGraphFlavor::kVg;

  // Branching graph: 1 -> 2, 1 -> 3
  //                       2 -> 4
  //                       3 -> 4
  imported.add_node({"1", "AAAA"});
  imported.add_node({"2", "CCCC"});
  imported.add_node({"3", "GGGG"});
  imported.add_node({"4", "TTTT"});

  imported.add_edge({"1", "2", false, false, "0", std::nullopt});
  imported.add_edge({"1", "3", false, false, "0", std::nullopt});
  imported.add_edge({"2", "4", false, false, "0", std::nullopt});
  imported.add_edge({"3", "4", false, false, "0", std::nullopt});

  // Path 1: 1 -> 2 -> 4
  ImportedPath path1;
  path1.name = "path1";
  path1.steps.push_back({"1", false});
  path1.steps.push_back({"2", false});
  path1.steps.push_back({"4", false});
  imported.add_path(path1);

  // Path 2: 1 -> 3 -> 4
  ImportedPath path2;
  path2.name = "path2";
  path2.steps.push_back({"1", false});
  path2.steps.push_back({"3", false});
  path2.steps.push_back({"4", false});
  imported.add_path(path2);

  PathGuidedTransform transform;
  AlnGraph result = transform.apply(imported, 0, 5);

  CHECK(result.pathCount() == 2);

  // Node 1 has two different successors (2 vs 3), so should be duplicated
  // Check that we have at least 2 variants of node 1 with different contexts
  int node1_variants = 0;
  for (std::size_t i = 0; i < result.nodeCount(); ++i) {
    const AlnNode& node = result.node(i);
    if (node.original_id == "1" && !node.is_reverse) {
      node1_variants++;
    }
  }

  CHECK(node1_variants >= 2);  // Should have duplicated node 1

  auto stats = transform.getStats();
  CHECK(stats.uncovered_node_count == 0);
}

//------------------------------------------------------------------------------
// Test 4: Graph with no embedded paths → all nodes uncovered
//------------------------------------------------------------------------------
TEST_CASE("PathGuidedTransform: no paths, all nodes uncovered") {
  ImportedGraph imported;
  imported.flavor = ImportedGraphFlavor::kVg;

  // Create nodes but NO paths
  imported.add_node({"1", "AAAA"});
  imported.add_node({"2", "CCCC"});
  imported.add_node({"3", "GGGG"});

  imported.add_edge({"1", "2", false, false, "0", std::nullopt});
  imported.add_edge({"2", "3", false, false, "0", std::nullopt});

  PathGuidedTransform transform;
  AlnGraph result = transform.apply(imported, 0, 5);

  auto stats = transform.getStats();
  CHECK(stats.original_node_count == 3);
  CHECK(stats.uncovered_node_count == 3);  // All nodes uncovered
  CHECK(result.pathCount() == 0);          // No paths created
}

//------------------------------------------------------------------------------
// Test 5: Path coverage statistics are accurate
//------------------------------------------------------------------------------
TEST_CASE("PathGuidedTransform: accurate coverage statistics") {
  ImportedGraph imported;
  imported.flavor = ImportedGraphFlavor::kVg;

  // 5 nodes, path covers only 3 of them
  imported.add_node({"1", "AAAA"});
  imported.add_node({"2", "CCCC"});
  imported.add_node({"3", "GGGG"});
  imported.add_node({"4", "TTTT"});  // Uncovered
  imported.add_node({"5", "NNNN"});  // Uncovered

  imported.add_edge({"1", "2", false, false, "0", std::nullopt});
  imported.add_edge({"2", "3", false, false, "0", std::nullopt});

  // Path covers only nodes 1, 2, 3
  ImportedPath path;
  path.name = "path1";
  path.steps.push_back({"1", false});
  path.steps.push_back({"2", false});
  path.steps.push_back({"3", false});
  imported.add_path(path);

  PathGuidedTransform transform;
  AlnGraph result = transform.apply(imported, 0, 5);

  auto stats = transform.getStats();
  CHECK(stats.original_node_count == 5);
  CHECK(stats.uncovered_node_count == 2);  // Nodes 4 and 5
  CHECK(result.pathCount() == 1);
}
