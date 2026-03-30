/**
 * bucket_indexer.hpp
 *
 * Bucket-partitioned indexer with concurrent producer-consumer.
 * Processes each node's interior exactly once (node-first), scatters
 * seeds into hash-partitioned buckets via lock-free queues, then
 * fills boundary gaps with path walks. Builds FlatSeedStore directly
 * from sorted buckets -- no HashSeedStore intermediate.
 *
 * Solves OOM for large genomes (hg38) by avoiding unordered_map
 * overhead (~80B/hash) and thread-store merge doubling.
 *
 * Related:
 *  - bucket_indexer.cpp
 *  - node_first_indexer.hpp (predecessor, uses HashSeedStore)
 *  - seed_store.hpp (FlatSeedStore output format)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/flat_graph.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/models/model.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"
#include "signal/seed_extractors/seed_extractor.hpp"

namespace piru::index {

struct BucketIndexConfig {
    /* Seed extraction */
    std::size_t seed_k{6};
    std::size_t seed_stride{1};
    double seed_freq_cutoff{1.0};    // keep bottom N% by frequency (1.0 = no filter)
    double seed_freq_cap{0.25};

    /* Bucket partitioning */

    // Number of hash buckets. 0 = auto (executor->max_concurrency(), or 1 if sequential).
    std::size_t num_buckets{0};

    // Backpressure: max total queued entries across all buckets before
    // producers yield. 16M entries * 16 bytes = 256 MB.
    std::size_t max_queued_entries{16 * 1024 * 1024};

    /* Parallelization */
    concurrency::Executor* executor{nullptr};
};

struct BucketIndexResult {
    // FlatSeedStore built directly from sorted buckets (no HashSeedStore)
    std::unique_ptr<FlatSeedStore> seed_store;

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
 * Bucket-partitioned indexing with concurrent producer-consumer.
 *
 * Phase 1 -- Node interiors:
 *   Producer threads (TBB parallel_for over nodes) squigglize using
 *   fast 2-bit rolling kmer + flat pore model array, quantize, extract
 *   seeds, scatter to ConcurrentQueue[hash % num_buckets].
 *   Consumer threads drain queues into flat vectors.
 *
 * Phase 2 -- Boundary fill:
 *   Walk paths, build k-mer context windows across node junctions,
 *   squigglize boundary, extract seeds, scatter to same queues.
 *   Also builds linearization coordinates.
 *
 * Phase 3 -- Finalize:
 *   Sort all entries by hash, build CSR arrays, apply global frequency
 *   filter. Construct FlatSeedStore directly.
 *
 * Signal processing matches path-walk behavior: raw pore model values
 * go directly to the fuzzy quantizer (no global normalization).
 */
BucketIndexResult bucketIndex(const FlatGraph& graph, const io::KmerModel& model,
                              const signal::FuzzyQuantizer& fuzzy_quantizer,
                              const signal::SeedExtractor& extractor,
                              const BucketIndexConfig& config = {});

}  // namespace piru::index
