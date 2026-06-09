/**
 * bucket_indexer.hpp
 *
 * Path-walk bucket indexer: walks embedded paths to extract seeds,
 * scatters to hash-partitioned buckets, finalizes each bucket
 * independently. Correct for both k-mer and minimizer extraction.
 *
 * Related:
 *  - bucket_indexer.cpp
 *  - bucket_seed_store.hpp (final index structure)
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
#include "signal/io/models/model.hpp"
#include "signal/tokenizers/tokenizer.hpp"
#include "signal/seed_extractors/seed_extractor.hpp"

namespace panomap::index {

struct BucketIndexConfig {
  /* Seed extraction params (must match extractor config) */
  std::size_t seed_k;  // Set by IndexPipelineConfig, no default.
  /* Signal processing */
  float diff_filter{0.0f};  // Diff filter threshold (0 = disabled)

  /* Parallelization */
  concurrency::Executor* executor{nullptr};
};

struct BucketIndexResult {
  // BucketSeedStore: per-bucket sorted keys + entries
  std::unique_ptr<SeedStore> seed_store;

  // Linear coordinates per node (indexed by node_id)
  std::vector<std::vector<LinearCoordinate>> linearization_coords;

  // Path lengths in base pairs (indexed by path_id)
  std::vector<std::size_t> path_lengths;

  // Stats
  std::size_t seeds_interior{0};
  std::size_t seeds_boundary{0};
  std::size_t seeds_total{0};
};

/**
 * Path-walk bucket indexing.
 *
 * Phase 1 -- Path walks:
 *   TBB parallel_for over paths. Squigglize full concatenated path,
 *   quantize, extract seeds, map positions to (node_id, offset),
 *   push to thread-local bucket buffers. Merge after all paths.
 *   Also builds linearization coordinates.
 *
 * Phase 2 -- Per-bucket finalization (parallel):
 *   Sort, dedup exact (hash, node_id, offset) duplicates from shared
 *   nodes, apply global frequency filter, build BucketSeedStore.
 *
 * .pirx serialization currently flattens BucketSeedStore to CSR
 * for compatibility. Bucket-native on-disk format is future work.
 */
BucketIndexResult bucketIndex(const FlatGraph& graph, const io::KmerModel& model,
                              const signal::Tokenizer& tokenizer,
                              const signal::SeedExtractor& extractor,
                              const BucketIndexConfig& config = {});

}  // namespace panomap::index
