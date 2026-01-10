// SPDX-License-Identifier: MIT
/**
 * Node-first indexing for the simple pipeline.
 *
 * Processes each node's interior exactly once with global normalization,
 * then fills boundary gaps during path walks. Faster than path-walk indexing
 * for graphs with high path sharing (conserved regions).
 *
 * Two-pass approach:
 *   Pass 1: Process node interiors (squigglize, accumulate global stats,
 *           normalize, fuzzy quantize, index positions with complete context)
 *   Pass 2: Walk paths to fill boundary gaps using getHashWindow()
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "index/aln_graph.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/models/model.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"
#include "signal/seed_extractors/seed_extractor.hpp"

namespace piru::index {

struct NodeFirstIndexConfig {
    // Seed extraction parameters
    std::size_t seed_k{6};      // Hash window size (fuzzy samples per seed)
    std::size_t seed_stride{1};
    double seed_filter{0.5};    // keep_least_frequent_fraction
};

struct NodeFirstIndexResult {
    // Seed store populated with (node_id, offset) entries
    std::unique_ptr<HashSeedStore> seed_store;

    // Linear coordinates per node (indexed by node_id)
    std::vector<std::vector<LinearCoordinate>> linearization_coords;

    // Path lengths in signal samples (indexed by path_id)
    std::vector<std::size_t> path_lengths;

    // Global normalization parameters (from pass 1)
    float global_mean{0.0f};
    float global_std{1.0f};

    // Stats
    std::size_t seeds_interior{0};   // Seeds from node interiors (pass 1)
    std::size_t seeds_boundary{0};   // Seeds from boundary fill (pass 2)
    std::size_t seeds_unique{0};     // After dedup
};

/**
 * Node-first indexing.
 *
 * Pass 1 - Node Interiors:
 *   For each node with length L (pore k-mer size k, hash window w):
 *   - Squigglize positions [0, L-k] -> L-k+1 raw k-mer values
 *   - Accumulate sum/count for global mean/std
 *   - After all nodes: compute global normalization params
 *   - Normalize, fuzzy quantize, index positions [0, L-k-w+1]
 *
 * Pass 2 - Boundary Fill:
 *   For each path, walk nodes and fill gaps:
 *   - Start at first unindexed position (L-k-w+2, or 0 if node too small)
 *   - For each position, call getHashWindow() to get fuzzy samples
 *   - Hash and index until end of node
 *   - Record linearization coordinates
 *
 * @param graph AlnGraph from simpleExpand() (2x nodes: forward + reverse)
 * @param model Pore model for squigglization (k-mer -> expected signal)
 * @param fuzzy_quantizer Converts normalized signal to discrete tokens
 * @param extractor Extracts seeds from fuzzy-quantized tokens
 * @param config Seed extraction parameters
 *
 * @return NodeFirstIndexResult with seed_store, linearization_coords,
 *         global normalization params, and stats
 */
NodeFirstIndexResult nodeFirstIndex(
    const AlnGraph& graph,
    const io::KmerModel& model,
    const signal::FuzzyQuantizer& fuzzy_quantizer,
    const signal::SeedExtractor& extractor,
    const NodeFirstIndexConfig& config = {});

}  // namespace piru::index
