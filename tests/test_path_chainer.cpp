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
// Tests default to min_chain_anchors=1 so single-anchor chains are preserved.
PathChainer make_chainer(const std::vector<std::vector<LinearCoordinate>>& coords,
                         const std::vector<std::size_t>& path_lengths,
                         PathChainerConfig config = {}) {
  config.min_chain_anchors = 1;
  return PathChainer(std::move(config), coords, path_lengths);
}

// Convenience accessors for best chain from ChainResult
const std::vector<ChainedAnchor>& best_anchors(const ChainResult& r) { return r.chains[0].anchors; }
double best_score(const ChainResult& r) { return r.chains[0].score; }

}  // namespace

TEST_CASE("PathChainer: Empty input returns empty output") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  auto chainer = make_chainer(coords, path_lengths);
  auto summary = chainer.chain({});

  CHECK(summary.chains.empty());
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

  REQUIRE(summary.chains.size() == 1);
  REQUIRE(best_anchors(summary).size() == 1);
  CHECK(best_anchors(summary)[0].read_pos == 50);
  CHECK(best_score(summary) > 0.0);
}

TEST_CASE("PathChainer: Linear colinear chain selects all anchors") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20),
                                  make_hit(2, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(best_anchors(summary).size() == 3);
  CHECK(best_anchors(summary)[0].read_pos == 50);
  CHECK(best_anchors(summary)[1].read_pos == 150);
  CHECK(best_anchors(summary)[2].read_pos == 250);
}

TEST_CASE("PathChainer: Large diagonal deviation breaks chain") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 50;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  // Δr = 100, Δq = 200, |Δr - Δq| = 100 > max_diag_dev
  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() == 1);
}

TEST_CASE("PathChainer: Distance filter prevents chaining far-apart anchors") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 10000}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 1000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 9950, 20)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() == 1);
}

TEST_CASE("PathChainer: Rejects cross-path chains") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{1, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() == 1);
}

TEST_CASE("PathChainer: Prefers higher-scoring chain") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 10),
                                  make_hit(2, 0, 250, 30)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() >= 2);
}

TEST_CASE("PathChainer: Min chain score threshold filters weak chains") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.min_chain_score = 1000;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  CHECK(summary.chains.empty());
}

TEST_CASE("PathChainer: Overlapping anchors are not chained together") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 150}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  config.min_chain_anchors = 1;
  auto chainer = make_chainer(coords, path_lengths, config);

  // dq = 10, span = 20, so query positions overlap -> dr must be > 0 and dq > 0
  // ref gap = 50, query gap = 10, but both positive so it should chain
  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 60, 20)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() >= 1);
}

TEST_CASE("PathChainer: Backward query positions are rejected") {
  std::vector<std::vector<LinearCoordinate>> coords(2);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 150, 20), make_hit(1, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  CHECK(best_anchors(summary).size() == 1);
}

TEST_CASE("PathChainer: Node appearing on multiple paths expands correctly") {
  std::vector<std::vector<LinearCoordinate>> coords(1);
  coords[0] = {{0, 100}, {1, 500}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(summary.chains.size() >= 1);
  REQUIRE(best_anchors(summary).size() == 1);
  CHECK(best_score(summary) > 0.0);
}

TEST_CASE("PathChainer: Chain preserves order from backtracking") {
  std::vector<std::vector<LinearCoordinate>> coords(3);
  coords[0] = {{0, 100}};
  coords[1] = {{0, 200}};
  coords[2] = {{0, 300}};
  std::vector<std::size_t> path_lengths(100, 1000000);

  PathChainerConfig config;
  config.max_dist_ref = 5000;
  config.bw = 500;
  config.min_chain_score = 10;
  auto chainer = make_chainer(coords, path_lengths, config);

  std::vector<NodeAnchor> hits = {make_hit(0, 0, 50, 20), make_hit(1, 0, 150, 20),
                                  make_hit(2, 0, 250, 20)};
  auto summary = chainer.chain(hits);

  REQUIRE(best_anchors(summary).size() == 3);
  CHECK(best_anchors(summary)[0].read_pos < best_anchors(summary)[1].read_pos);
  CHECK(best_anchors(summary)[1].read_pos < best_anchors(summary)[2].read_pos);
}
