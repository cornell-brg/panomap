// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "mapping/seed_merger.hpp"

using namespace piru::mapping;

namespace {

// Helper to create a NodeAnchor for testing
NodeAnchor make_hit(std::size_t node_id, std::size_t offset, std::size_t query_pos,
                    std::size_t span) {
  NodeAnchor hit;
  hit.node_id = static_cast<std::uint32_t>(node_id);
  hit.offset = static_cast<std::uint32_t>(offset);
  hit.read_pos = static_cast<std::uint32_t>(query_pos);
  hit.span = static_cast<std::uint16_t>(span);
  hit.length = static_cast<std::uint16_t>(span);
  return hit;
}

}  // namespace

TEST_CASE("SeedMerger: Empty input returns empty output") {
  SeedMergerConfig config{.merge_tolerance = 0};
  auto merged = SeedMerger::merge({}, config);
  CHECK(merged.empty());
}

TEST_CASE("SeedMerger: Single hit returns unchanged") {
  SeedMergerConfig config{.merge_tolerance = 0};
  std::vector<NodeAnchor> hits = {make_hit(1, 10, 20, 15)};

  auto merged = SeedMerger::merge(hits, config);

  REQUIRE(merged.size() == 1);
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].offset == 10);
  CHECK(merged[0].read_pos == 20);
  CHECK(merged[0].span == 15);
}

TEST_CASE("SeedMerger: Gap with tolerance=0 keeps seeds separate") {
  SeedMergerConfig config{.merge_tolerance = 0};

  // Two seeds on same node with gap=5 in both query and reference
  // Seed 1: query 100-120, ref 10-30
  // Seed 2: query 125-145, ref 35-55 (gap=5 in both)
  std::vector<NodeAnchor> hits = {
      make_hit(1, 10, 100, 20),  // Ends at query=120, ref=30
      make_hit(1, 35, 125, 20)   // Starts at query=125, ref=35 (gap=5 > tolerance=0)
  };

  auto merged = SeedMerger::merge(hits, config);

  // Should NOT merge with tolerance=0 (gap=5 > 0)
  REQUIRE(merged.size() == 2);
}

TEST_CASE("SeedMerger: Perfect overlap merges with tolerance=0") {
  SeedMergerConfig config{.merge_tolerance = 0};

  // Two identical seeds (perfect overlap)
  std::vector<NodeAnchor> hits = {make_hit(1, 10, 100, 20), make_hit(1, 10, 100, 20)};

  auto merged = SeedMerger::merge(hits, config);

  // Should merge (gap = 0)
  REQUIRE(merged.size() == 1);
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].offset == 10);
  CHECK(merged[0].read_pos == 100);
  CHECK(merged[0].span == 20);  // Unchanged since they overlap perfectly
}

TEST_CASE("SeedMerger: Gap within tolerance merges seeds") {
  SeedMergerConfig config{.merge_tolerance = 10};

  // Two seeds with gap=5 in both query and reference
  std::vector<NodeAnchor> hits = {
      make_hit(1, 10, 100, 20),  // Ends at query=120, ref=30
      make_hit(1, 15, 105, 20)   // Starts at query=105, ref=15 (gap=5 in both)
  };

  auto merged = SeedMerger::merge(hits, config);

  // Should merge (gap=5 <= tolerance=10)
  REQUIRE(merged.size() == 1);
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].offset == 10);     // Use earlier offset
  CHECK(merged[0].read_pos == 100);  // Use earlier query pos
  CHECK(merged[0].span == 25);       // Covers 100 to 125 (105+20)
}

TEST_CASE("SeedMerger: Gap exceeds tolerance keeps seeds separate") {
  SeedMergerConfig config{.merge_tolerance = 5};

  // Two seeds with gap=10 in query (exceeds tolerance)
  // Seed 1: query 100-120, ref 10-30
  // Seed 2: query 130-150, ref 40-60 (gap=10 in query, ref)
  std::vector<NodeAnchor> hits = {
      make_hit(1, 10, 100, 20),  // Ends at query=120, ref=30
      make_hit(1, 40, 130, 20)   // Starts at query=130, ref=40 (gap=10 > tolerance=5)
  };

  auto merged = SeedMerger::merge(hits, config);

  // Should NOT merge (gap=10 > tolerance=5)
  REQUIRE(merged.size() == 2);
}

TEST_CASE("SeedMerger: Different nodes do not merge") {
  SeedMergerConfig config{.merge_tolerance = 100};  // Large tolerance

  // Two seeds on different nodes (but close in query positions)
  std::vector<NodeAnchor> hits = {
      make_hit(1, 10, 100, 20), make_hit(2, 10, 105, 20)  // Different node_id
  };

  auto merged = SeedMerger::merge(hits, config);

  // Should NOT merge (different nodes)
  REQUIRE(merged.size() == 2);
  CHECK(merged[0].node_id == 1);
  CHECK(merged[1].node_id == 2);
}

TEST_CASE("SeedMerger: Merges multiple adjacent seeds (A+B+C)") {
  SeedMergerConfig config{.merge_tolerance = 5};

  // Three seeds that can all be merged
  std::vector<NodeAnchor> hits = {
      make_hit(1, 10, 100, 10),  // Ends at query=110, ref=20
      make_hit(1, 12, 105, 10),  // Gap=5, ends at query=115, ref=22
      make_hit(1, 14, 108, 10)   // Gap=3 (from previous), ends at query=118, ref=24
  };

  auto merged = SeedMerger::merge(hits, config);

  // All three should merge into one
  REQUIRE(merged.size() == 1);
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].offset == 10);     // Earliest offset
  CHECK(merged[0].read_pos == 100);  // Earliest query pos
  CHECK(merged[0].span == 18);       // Covers 100 to 118
}

TEST_CASE("SeedMerger: Handles unsorted input correctly") {
  SeedMergerConfig config{.merge_tolerance = 5};

  // Provide hits in unsorted order
  std::vector<NodeAnchor> hits = {
      make_hit(2, 20, 200, 10), make_hit(1, 10, 100, 10),
      make_hit(1, 12, 105, 10)  // Should merge with previous seed on node 1
  };

  auto merged = SeedMerger::merge(hits, config);

  // Should produce: one merged hit on node 1, one hit on node 2
  REQUIRE(merged.size() == 2);

  // First result should be merged node 1 seeds (sorted by node_id)
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].read_pos == 100);
  CHECK(merged[0].span == 15);  // Covers 100 to 115

  // Second result should be node 2 seed
  CHECK(merged[1].node_id == 2);
  CHECK(merged[1].read_pos == 200);
}

TEST_CASE("SeedMerger: Overlapping seeds merge correctly") {
  SeedMergerConfig config{.merge_tolerance = 20};

  // Two seeds with significant overlap
  // Seed 1: query 100-120 (span=20)
  // Seed 2: query 110-130 (span=20, overlap=10)
  std::vector<NodeAnchor> hits = {make_hit(1, 10, 100, 20), make_hit(1, 15, 110, 20)};

  auto merged = SeedMerger::merge(hits, config);

  // Should merge into one seed covering 100-130
  REQUIRE(merged.size() == 1);
  CHECK(merged[0].read_pos == 100);
  CHECK(merged[0].span == 30);  // 100 to 130
}

TEST_CASE("SeedMerger: Multiple clusters remain separate") {
  SeedMergerConfig config{.merge_tolerance = 5};

  // Two clusters of mergeable seeds with large gap between clusters
  std::vector<NodeAnchor> hits = {// Cluster 1 on node 1
                                  make_hit(1, 10, 100, 10), make_hit(1, 12, 105, 10),

                                  // Cluster 2 on node 1 (far away in query space)
                                  make_hit(1, 50, 500, 10), make_hit(1, 52, 505, 10),

                                  // Cluster 3 on node 2
                                  make_hit(2, 10, 100, 10), make_hit(2, 12, 105, 10)};

  auto merged = SeedMerger::merge(hits, config);

  // Should produce 3 merged hits (one per cluster)
  REQUIRE(merged.size() == 3);

  // Check cluster 1 (node 1, query ~100)
  CHECK(merged[0].node_id == 1);
  CHECK(merged[0].read_pos == 100);
  CHECK(merged[0].span == 15);

  // Check cluster 2 (node 1, query ~500)
  CHECK(merged[1].node_id == 1);
  CHECK(merged[1].read_pos == 500);
  CHECK(merged[1].span == 15);

  // Check cluster 3 (node 2, query ~100)
  CHECK(merged[2].node_id == 2);
  CHECK(merged[2].read_pos == 100);
  CHECK(merged[2].span == 15);
}
