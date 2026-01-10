// SPDX-License-Identifier: MIT

#include "mapping/dp_chain_clusterer.hpp"

#include <doctest/doctest.h>

#include "mapping/anchor_expander.hpp"

using namespace piru::mapping;
using namespace piru::index;

namespace {

// Helper to create a SeedHitRecord for testing
SeedHitRecord make_hit(std::size_t node_id, std::size_t offset, std::size_t query_pos,
                       std::size_t span) {
    SeedHitRecord hit;
    hit.target.node_id = node_id;
    hit.target.offset = offset;
    hit.target.length = span;
    hit.read_pos = query_pos;
    hit.span = span;
    hit.hash = 0;
    return hit;
}

}  // namespace

TEST_CASE("DPChainClusterer: Empty input returns empty output") {
    // Setup linearization (1 node, 1 path)
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    auto anchors = expander.expand({});

    DPChainClustererConfig config;
    DPChainClusterer clusterer(config);

    auto summary = clusterer.cluster(anchors);

    CHECK(summary.anchors.empty());
    CHECK(summary.score == 0.0);
}

TEST_CASE("DPChainClusterer: Single seed produces single anchor chain") {
    // Node 0 on path 0 at position 100
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.min_chain_score = 10;  // Lower threshold

    DPChainClusterer clusterer(config);

    // Single seed hit
    std::vector<SeedHitRecord> hits = {make_hit(0, 0, 50, 20)};
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    REQUIRE(summary.anchors.size() == 1);
    CHECK(summary.anchors[0].read_pos == 50);
    CHECK(summary.score > 0.0);
}

TEST_CASE("DPChainClusterer: Linear colinear chain selects all anchors") {
    // Three nodes in a line on same path
    std::vector<std::vector<LinearCoordinate>> coords(3);
    coords[0] = {{0, 100}};   // Node 0 at ref 100
    coords[1] = {{0, 200}};   // Node 1 at ref 200
    coords[2] = {{0, 300}};   // Node 2 at ref 300

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Three colinear seeds (diagonal: ref +100, query +100)
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),    // ref 100, query 50
        make_hit(1, 0, 150, 20),   // ref 200, query 150 (Δr=100, Δq=100, diag=0)
        make_hit(2, 0, 250, 20)    // ref 300, query 250 (Δr=100, Δq=100, diag=0)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // All three should be chained
    REQUIRE(summary.anchors.size() == 3);
    CHECK(summary.anchors[0].read_pos == 50);
    CHECK(summary.anchors[1].read_pos == 150);
    CHECK(summary.anchors[2].read_pos == 250);
}

TEST_CASE("DPChainClusterer: Large diagonal deviation breaks chain") {
    // Two nodes on same path
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 200}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 50;  // Strict diagonal constraint
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Two seeds with large diagonal deviation
    // Δr = 100, Δq = 200, |Δr - Δq| = 100 > max_diag_dev
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),    // ref 100, query 50
        make_hit(1, 0, 250, 20)    // ref 200, query 250 (diagonal dev = 100)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should only select one anchor (cannot chain due to diagonal constraint)
    CHECK(summary.anchors.size() == 1);
}

TEST_CASE("DPChainClusterer: Distance filter prevents chaining far-apart anchors") {
    // Two nodes on same path
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 10000}};  // Far away

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 1000;  // Strict distance limit
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Two seeds far apart (Δr = 9900 > max_dist)
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),
        make_hit(1, 0, 9950, 20)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should only select one anchor (too far to chain)
    CHECK(summary.anchors.size() == 1);
}

TEST_CASE("DPChainClusterer: Rejects cross-path chains") {
    // Node 0 on path 0, Node 1 on path 1
    // Cross-path chaining not supported (see DEV033 for planned alias-based approach)
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{1, 200}};  // Different path

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Two seeds on different paths
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),
        make_hit(1, 0, 150, 20)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should only select one anchor (cannot cross paths)
    CHECK(summary.anchors.size() == 1);
}

TEST_CASE("DPChainClusterer: Prefers higher-scoring chain") {
    // Three nodes on same path
    std::vector<std::vector<LinearCoordinate>> coords(3);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 200}};
    coords[2] = {{0, 300}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;
    config.anchor_weight = 1.0;

    DPChainClusterer clusterer(config);

    // Three seeds: 0→1 is good, 0→2 skips 1 but has larger span
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),    // Score = 20
        make_hit(1, 0, 150, 10),   // Score = 10 (smaller)
        make_hit(2, 0, 250, 30)    // Score = 30 (larger)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should select chain that maximizes total score
    // Best chain is likely all three: 20 + 10 + 30 - penalties
    CHECK(summary.anchors.size() >= 2);
}

TEST_CASE("DPChainClusterer: Min chain score threshold filters weak chains") {
    // Single node
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.min_chain_score = 1000;  // Very high threshold
    config.anchor_weight = 1.0;

    DPChainClusterer clusterer(config);

    // Single small seed (score = 20 < threshold)
    std::vector<SeedHitRecord> hits = {make_hit(0, 0, 50, 20)};
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should return empty (chain score below threshold)
    CHECK(summary.anchors.empty());
}

TEST_CASE("DPChainClusterer: Overlapping anchors incur penalty") {
    // Two nodes on same path
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 150}};  // Close together (potential overlap)

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.overlap_penalty_factor = 2.0;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Two seeds that overlap in reference space
    // Seed 0: ref 100-120, query 50-70
    // Seed 1: ref 150-170, query 60-80 (overlap in query)
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),
        make_hit(1, 0, 60, 20)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should still chain (overlap is allowed, just penalized)
    CHECK(summary.anchors.size() >= 1);
}

TEST_CASE("DPChainClusterer: Backward query positions are rejected") {
    // Two nodes on same path
    std::vector<std::vector<LinearCoordinate>> coords(2);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 200}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Second seed has earlier query position (backward)
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 150, 20),   // query 150
        make_hit(1, 0, 50, 20)     // query 50 (backward!)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should only select one (cannot chain backward in query)
    CHECK(summary.anchors.size() == 1);
}

TEST_CASE("DPChainClusterer: Node appearing on multiple paths expands correctly") {
    // Node 0 appears on two paths
    // Cross-path chaining not supported, so only one path's anchor is selected
    std::vector<std::vector<LinearCoordinate>> coords(1);
    coords[0] = {{0, 100}, {1, 500}};  // Same node on path 0 and path 1

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Single seed hit (expands to 2 anchors on different paths)
    std::vector<SeedHitRecord> hits = {make_hit(0, 0, 50, 20)};
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    // Should return 1 anchor (from best path)
    REQUIRE(summary.anchors.size() == 1);
    CHECK(summary.score > 0.0);
}

TEST_CASE("DPChainClusterer: Chain preserves order from backtracking") {
    // Three nodes in order
    std::vector<std::vector<LinearCoordinate>> coords(3);
    coords[0] = {{0, 100}};
    coords[1] = {{0, 200}};
    coords[2] = {{0, 300}};

    std::vector<std::size_t> path_lengths(100, 1000000);  // Large lengths to not filter any anchors
    PathWalkExpander expander(coords, path_lengths);
    DPChainClustererConfig config;
    config.max_dist = 5000;
    config.max_diag_dev = 500;
    config.min_chain_score = 10;

    DPChainClusterer clusterer(config);

    // Seeds in order
    std::vector<SeedHitRecord> hits = {
        make_hit(0, 0, 50, 20),
        make_hit(1, 0, 150, 20),
        make_hit(2, 0, 250, 20)
    };
    auto anchors = expander.expand(hits);

    auto summary = clusterer.cluster(anchors);

    REQUIRE(summary.anchors.size() == 3);
    // Check anchors are in forward query order
    CHECK(summary.anchors[0].read_pos < summary.anchors[1].read_pos);
    CHECK(summary.anchors[1].read_pos < summary.anchors[2].read_pos);
}
