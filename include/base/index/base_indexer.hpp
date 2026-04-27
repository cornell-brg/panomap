/**
 * base_indexer.hpp
 *
 * Path-walk minimizer indexer for base mode. Walks GFA paths, extracts
 * directional 2-bit k-mers (k, w configurable), hashes with minimap2's
 * hash64 mixer, scatters minimizer-selected hits into hash-partitioned
 * buckets, finalizes into a BucketSeedStore.
 *
 * Mirrors the structure of signal/index/bucket_indexer.cpp but without
 * the squigglize/tokenize/extract pipeline -- minimizers are computed
 * directly from the base sequences exposed by FlatGraph.
 *
 * Related:
 *  - base_indexer.cpp
 *  - core/index/bucket_seed_store.hpp (final index structure)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "core/concurrency/executor.hpp"
#include "core/index/flat_graph.hpp"
#include "core/index/linearizer.hpp"
#include "core/index/seed_store.hpp"

namespace piru::base {

struct BaseBucketIndexConfig {
  std::size_t k{15};       // k-mer length
  std::size_t w{10};       // minimizer window
  concurrency::Executor* executor{nullptr};
};

struct BaseBucketIndexResult {
  std::unique_ptr<piru::index::SeedStore> seed_store;
  std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
  std::vector<std::size_t> path_lengths;
  std::size_t seeds_total{0};
};

// Build a base-mode minimizer index from a directional FlatGraph. Walks
// each path's bases, computes hash64 over rolling k-mers (skipping any
// k-mer containing an N base), sliding-window minimizer with leftmost
// tiebreak. Hits are mapped back to (node_id, local_offset) via path
// boundaries before bucket scatter.
BaseBucketIndexResult bucketIndexBase(const piru::index::FlatGraph& graph,
                                      const BaseBucketIndexConfig& config = {});

}  // namespace piru::base
