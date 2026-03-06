// SPDX-License-Identifier: MIT
/**
 * Unified path-walk indexing for the simple pipeline.
 *
 * Combines squigglization, linearization, and seed extraction into a single
 * path-walking pass (conceptually - two passes needed for normalization stats).
 *
 * This replaces the separate squigglize -> linearize -> seed_builder stages
 * of the classic pipeline.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/aln_graph.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/models/model.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"
#include "signal/seed_extractors/seed_extractor.hpp"

namespace piru::index {

struct PathWalkIndexConfig {
  // Seed extraction parameters
  std::size_t seed_k{6};
  std::size_t seed_stride{1};
  double seed_freq_cutoff{0.5};  // threshold percentile (above this -> subsample)
  double seed_freq_cap{0.55};    // subsample cap percentile (target for subsampled seeds)

  // Debug: dump per-path normalization stats to file (empty = disabled)
  std::string dump_norm_stats_path;

  // Parallelization
  // If non-null, enables parallel indexing using this executor.
  // Caller owns the executor lifetime.
  concurrency::Executor* executor{nullptr};
};

struct PathWalkIndexResult {
  // Seed store populated with (node_id, offset) entries
  std::unique_ptr<HashSeedStore> seed_store;

  // Linear coordinates per node (indexed by node_id)
  // Each node may have multiple coordinates if it appears on multiple paths
  std::vector<std::vector<LinearCoordinate>> linearization_coords;

  // Path lengths in signal samples (indexed by path_id)
  std::vector<std::size_t> path_lengths;

  // Stats
  std::size_t total_path_length{0};  // Total signal samples across all paths
  std::size_t seeds_extracted{0};    // Before dedup
  std::size_t seeds_unique{0};       // After dedup
};

/**
 * Unified path-walk indexing.
 *
 * Walks all paths in the graph exactly twice:
 *   Pass 1: Squigglize each node's sequence, accumulate normalization stats
 *   Pass 2: Normalize, fuzzy quantize, extract seeds, track linear coordinates
 *
 * @param graph AlnGraph from simpleExpand() (2x nodes: forward + reverse)
 * @param model Pore model for squigglization (k-mer -> expected signal)
 * @param fuzzy_quantizer Converts normalized signal to discrete tokens
 * @param extractor Extracts seeds from fuzzy-quantized tokens
 * @param config Seed extraction parameters
 *
 * @return seed_store (hash -> [(node_id, offset), ...] with dedup) and
 *         linearization_coords (per-node linear coordinates from path walks)
 *
 * @note No SignalStore output - the simple pipeline doesn't store per-node signals.
 */
PathWalkIndexResult pathWalkIndex(const AlnGraph& graph, const io::KmerModel& model,
                                  const signal::FuzzyQuantizer& fuzzy_quantizer,
                                  const signal::SeedExtractor& extractor,
                                  const PathWalkIndexConfig& config = {});

}  // namespace piru::index
