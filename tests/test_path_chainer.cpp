// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "mapping/path_chainer.hpp"

using namespace piru::mapping;
using namespace piru::index;

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

// Helper to create a PathChainer with linearization data.
PathChainer make_chainer(const std::vector<std::vector<LinearCoordinate>>& coords,
                       const std::vector<std::size_t>& path_lengths, PathChainerConfig config = {}) {
  return PathChainer(std::move(config), coords, path_lengths);
}

}  // namespace

TEST_CASE("PathChainer: Empty input returns empty output") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  auto chainer = make_chainer(coords, path_lengths);
  auto summary = chainer.chain({});

  CHECK(summary.anchors.empty());
  CHECK(summary.score == 0.0);
}

TEST_CASE("PathChainer: Single seed produces single anchor chain") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(summary.anchors.size() == 1);
  CHECK(summary.anchors[0].read_pos == 50);
  CHECK(summary.score > 0.0);
}

TEST_CASE("PathChainer: Linear colinear chain selects all anchors") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20),
                                  make_hit(2, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(summary.anchors.size() == 3);
  CHECK(summary.anchors[0].read_pos == 50);
  CHECK(summary.anchors[1].read_pos == 150);
  CHECK(summary.anchors[2].read_pos == 250);
}

TEST_CASE("PathChainer: Large diagonal deviation breaks chain") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 50;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  // Δr = 100, Δq = 200, |Δr - Δq| = 100 > max_diag_dev
  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() == 1);
}

TEST_CASE("PathChainer: Distance filter prevents chaining far-apart anchors") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 10000}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 1000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 9950, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() == 1);
}

TEST_CASE("PathChainer: Rejects cross-path chains") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{1, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() == 1);
}

TEST_CASE("PathChainer: Prefers higher-scoring chain") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  config.anchor_weight = 1.0;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 10),
                                  make_hit(2, 0, 250, 30)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() >= 2);
}

TEST_CASE("PathChainer: Min chain score threshold filters weak chains") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.min_chain_score = 1000;
  config.anchor_weight = 1.0;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.empty());
}

TEST_CASE("PathChainer: Overlapping anchors incur penalty") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 150}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.overlap_penalty_factor = 2.0;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 60, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() >= 1);
}

TEST_CASE("PathChainer: Backward query positions are rejected") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 150, 20), make_hit(1, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.anchors.size() == 1);
}

TEST_CASE("PathChainer: Node appearing on multiple paths expands correctly") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}, {1, 500}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(summary.anchors.size() == 1);
  CHECK(summary.score > 0.0);
}

TEST_CASE("PathChainer: Chain preserves order from backtracking") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist = 5000;
  config.max_diag_dev = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20),
                                  make_hit(2, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(summary.anchors.size() == 3);
  CHECK(summary.anchors[0].read_pos < summary.anchors[1].read_pos);
  CHECK(summary.anchors[1].read_pos < summary.anchors[2].read_pos);
}
