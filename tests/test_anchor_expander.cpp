// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "mapping/anchor_expander.hpp"

using namespace piru::mapping;
using namespace piru::index;

namespace {

// Helper to create a NodeAnchor for testing
NodeAnchor make_hit(std::size_t node_id, std::size_t offset, std::size_t query_pos,
                       std::size_t span) {
    NodeAnchor hit;
    hit.target.node_id = node_id;
    hit.target.offset = offset;
    hit.target.length = span;
    hit.read_pos = query_pos;
    hit.span = span;
    hit.hash = 0;
    return hit;
}

}  // namespace

TEST_CASE("PathWalkExpander: Empty input returns empty output") {
    // Node 0 has one coordinate on path 0 at position 100
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    auto anchors = expander.expand({});

    CHECK(anchors.empty());
}

TEST_CASE("PathWalkExpander: Single seed on single path produces one anchor") {
    // Node 0 has one coordinate on path 0 at position 100
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed hit at node 0, offset 10, query position 50, span 20
    std::vector<NodeAnchor> hits = {make_hit(0, 10, 50, 20)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 1);
    CHECK(anchors[0].query_pos == 50);
    CHECK(anchors[0].ref_coord == 110);  // 100 + 10 (node start + offset)
    CHECK(anchors[0].length == 20);
    CHECK(anchors[0].path_id == 0);
    CHECK(anchors[0].node_id == 0);
    CHECK(anchors[0].node_offset == 10);
}

TEST_CASE("PathWalkExpander: Single seed on multiple paths produces multiple anchors") {
    // Node 1 appears on two paths:
    // - Path 0 at position 200
    // - Path 1 at position 500
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {};  // Node 0 has no coordinates
    coords[1] = {{0, 200}, {1, 500}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed hit at node 1, offset 5, query position 100, span 15
    std::vector<NodeAnchor> hits = {make_hit(1, 5, 100, 15)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 2);

    // First anchor on path 0
    CHECK(anchors[0].query_pos == 100);
    CHECK(anchors[0].ref_coord == 205);  // 200 + 5
    CHECK(anchors[0].length == 15);
    CHECK(anchors[0].path_id == 0);
    CHECK(anchors[0].node_id == 1);
    CHECK(anchors[0].node_offset == 5);

    // Second anchor on path 1
    CHECK(anchors[1].query_pos == 100);
    CHECK(anchors[1].ref_coord == 505);  // 500 + 5
    CHECK(anchors[1].length == 15);
    CHECK(anchors[1].path_id == 1);
    CHECK(anchors[1].node_id == 1);
    CHECK(anchors[1].node_offset == 5);
}

TEST_CASE("PathWalkExpander: Seed with no coordinates produces no anchors") {
    // Node 0 has coordinates, but node 1 does not
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {};  // Empty vector = no coordinates

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed hit at node 1 (no coordinates)
    std::vector<NodeAnchor> hits = {make_hit(1, 0, 50, 20)};

    auto anchors = expander.expand(hits);

    // Should produce no anchors (node has no linearization)
    CHECK(anchors.empty());
}

TEST_CASE("PathWalkExpander: Offset handling adds to ref_coord") {
    // Node 0 at position 1000 on path 0
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 1000}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed at offset 50
    std::vector<NodeAnchor> hits = {make_hit(0, 50, 100, 20)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 1);
    CHECK(anchors[0].ref_coord == 1050);  // 1000 + 50
}

TEST_CASE("PathWalkExpander: Zero offset produces ref_coord equal to node start") {
    // Node 0 at position 2000 on path 0
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 2000}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed at offset 0
    std::vector<NodeAnchor> hits = {make_hit(0, 0, 100, 20)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 1);
    CHECK(anchors[0].ref_coord == 2000);  // Exactly at node start
}

TEST_CASE("PathWalkExpander: Multiple seeds produce multiple anchors") {
    // Node 0 on path 0 at 100, Node 1 on path 0 at 200
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 200}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    std::vector<NodeAnchor> hits = {
        make_hit(0, 10, 50, 20),  // First seed on node 0
        make_hit(1, 5, 70, 15)    // Second seed on node 1
    };

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 2);

    // First anchor from node 0
    CHECK(anchors[0].query_pos == 50);
    CHECK(anchors[0].ref_coord == 110);
    CHECK(anchors[0].node_id == 0);

    // Second anchor from node 1
    CHECK(anchors[1].query_pos == 70);
    CHECK(anchors[1].ref_coord == 205);
    CHECK(anchors[1].node_id == 1);
}

TEST_CASE(
    "PathWalkExpander: Node appearing on same path multiple times produces multiple anchors") {
    // Node 0 appears twice on path 0 (e.g., in a cycle)
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}, {0, 500}};  // Same path_id, different positions

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Single seed hit at node 0
    std::vector<NodeAnchor> hits = {make_hit(0, 10, 50, 20)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 2);

    // First occurrence
    CHECK(anchors[0].query_pos == 50);
    CHECK(anchors[0].ref_coord == 110);  // 100 + 10
    CHECK(anchors[0].path_id == 0);

    // Second occurrence
    CHECK(anchors[1].query_pos == 50);
    CHECK(anchors[1].ref_coord == 510);  // 500 + 10
    CHECK(anchors[1].path_id == 0);
}

TEST_CASE("PathWalkExpander: Invalid node_id is skipped") {
    // Only node 0 has coordinates
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed hit at node 5 (out of bounds)
    std::vector<NodeAnchor> hits = {make_hit(5, 0, 50, 20)};

    auto anchors = expander.expand(hits);

    // Should skip invalid node
    CHECK(anchors.empty());
}

TEST_CASE("PathWalkExpander: Mix of valid and invalid nodes") {
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {};  // No coordinates

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    std::vector<NodeAnchor> hits = {
        make_hit(0, 10, 50, 20),  // Valid (has coords)
        make_hit(1, 5, 70, 15),   // Invalid (no coords)
        make_hit(5, 0, 90, 10)    // Invalid (out of bounds)
    };

    auto anchors = expander.expand(hits);

    // Only first seed should produce anchor
    REQUIRE(anchors.size() == 1);
    CHECK(anchors[0].query_pos == 50);
    CHECK(anchors[0].ref_coord == 110);
    CHECK(anchors[0].node_id == 0);
}

TEST_CASE("PathWalkExpander: Seed length propagates to anchor") {
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);

    // Seed with specific length (e.g., after merging)
    std::vector<NodeAnchor> hits = {make_hit(0, 0, 100, 50)};

    auto anchors = expander.expand(hits);

    REQUIRE(anchors.size() == 1);
    CHECK(anchors[0].length == 50);  // Length preserved
}
