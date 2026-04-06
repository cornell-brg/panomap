/**
 * index_pipeline.hpp
 *
 * Indexing pipeline: transforms an imported graph into a searchable
 * index (GraphStore + SeedStore + linearization coordinates).
 *
 * IndexPipelineConfig is the single source of truth for all indexing
 * parameters -- CLI defaults in index.cpp mirror these values.
 *
 * Related:
 *  - index_pipeline.cpp
 *  - node_first_indexer.cpp, path_walk_indexer.cpp  (backends)
 *  - serialization.hpp  (.pirx persistence)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "index/sort_1d.hpp"
#include "io/graphs/graph.hpp"
#include "io/models/model.hpp"

namespace piru::index {

/**
 * Configuration for the indexing pipeline.
 */
struct IndexPipelineConfig {
  /* Indexer backend */

  // Indexer backend. Only "bucket" is supported (path-walk generation +
  // bucket-finalized store). Legacy backends (node-first, path-walk) were
  // removed -- see git history if needed for reference.
  std::string indexer_backend{"bucket"};

  /* Signal processing */

  // Tokenizer: "rh2" (rawhash2)
  // Converts normalized signal to discrete tokens (4-bit = 16 values)
  std::string tokenizer{"rh2"};
  float tokenizer_fine_min{-2.0f};
  float tokenizer_fine_max{2.0f};
  float tokenizer_fine_range{0.4f};
  float tokenizer_diff{0.35f};  // Skip events within diff of last emitted (0 = disabled, RH2: 0.35)
  std::uint32_t tokenizer_n_bins{0};  // 0 = use 2^qbits = 16

  /* Seed extraction */

  std::string seed_type{"kmer"};    // "kmer" or "minimizer"
  std::size_t minimizer_window{5};  // minimizer window (w=1 -> rolling k-mer)
  std::size_t seed_k{8};            // tokens per seed hash
  std::size_t seed_stride{1};       // spacing between seeds (1 = dense)
  double seed_freq_cutoff{1.0};     // keep bottom N% by frequency (1.0 = no index-time filter)
  double seed_freq_cap{0.25};       // subsample cap for seeds above cutoff

  /* 1D sort (for SortChainer) */

  bool compute_1d_sort{true};   // compute 1D canonical coordinates (always on)
  Sort1DConfig sort_1d_config;  // SGD parameters

  /* Debug */

  std::string dump_norm_stats_path;  // dump per-path norm stats (path-walk only)

  /* Parallelization */

  // If nullptr, indexing runs sequentially.
  concurrency::Executor* executor{nullptr};

  // Note: qbits hardcoded to 4 in both tokenizer.hpp and seed extractor.
  // Do not change without coordinating both.
  //
  // Pore model k: determined by model (r9.4->k=6, r10.4->k=9), not configurable.
};

/**
 * Result of running the indexing pipeline (in-memory representation).
 */
struct IndexPipelineResult {
  std::unique_ptr<GraphStore> graph_store;
  std::unique_ptr<SeedStore> seed_store;
  std::vector<std::vector<LinearCoordinate>> linearization_coords;
  std::vector<float> node_1d_coords;  // 1D SGD positions for SortChainer (empty if not computed)
  std::vector<std::size_t> path_lengths;
  std::size_t pore_k{0};
  std::string model_name;
  std::string tokenizer;
};

/**
 * Run the full indexing pipeline on an imported graph.
 *
 * Stages: expand graph -> linearize -> squigglize -> quantize -> build seeds.
 */
// Takes imported by value so it can be freed after graph expansion.
IndexPipelineResult run_index_pipeline(io::ImportedGraph imported, const io::KmerModel& model,
                                       const IndexPipelineConfig& config);

}  // namespace piru::index
