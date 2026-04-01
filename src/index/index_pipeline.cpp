/**
 * index_pipeline.cpp
 *
 * Implements the indexing pipeline: graph expansion, signal generation,
 * quantization, and seed extraction via the bucket indexer.
 *
 * Legacy backends (node-first, path-walk) were removed. See git history
 * for their implementations if needed for reference.
 *
 * SPDX-License-Identifier: MIT
 */

#include "index/index_pipeline.hpp"

#include <chrono>
#include <memory>
#include <stdexcept>

#include "index/flat_graph.hpp"
#include "index/bucket_indexer.hpp"
#include "index/simple_expand.hpp"
#include "index/sort_1d.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "util/logging.hpp"

namespace piru::index {

IndexPipelineResult run_index_pipeline(io::ImportedGraph imported,
                                       const io::KmerModel& model,
                                       const IndexPipelineConfig& config) {
  auto stage_start = std::chrono::high_resolution_clock::now();

  /* 1. Expand graph (forward + reverse nodes) */

  auto flat_graph = simpleExpandFlat(imported);

  // Free imported graph ASCII sequences (~3.1 GB for hg38)
  imported = io::ImportedGraph{};

  auto stage_elapsed =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - stage_start)
          .count();
  LOG_INFO(
      "[1/2] Transformed graph to directional graph: " + std::to_string(flat_graph.nodeCount()) +
      " nodes, " + std::to_string(flat_graph.edgeCount()) + " edges, " +
      std::to_string(flat_graph.pathCount()) + " paths [" + std::to_string(stage_elapsed) + "s]");

  /* 2. Index (squigglize + linearize + seed extraction) */

  stage_start = std::chrono::high_resolution_clock::now();

  signal::FuzzyQuantizerConfig fuzzy_cfg;
  fuzzy_cfg.backend = config.fuzzy_quantizer;
  fuzzy_cfg.pore_model = model.name();
  fuzzy_cfg.fine_min = config.fuzzy_fine_min;
  fuzzy_cfg.fine_max = config.fuzzy_fine_max;
  fuzzy_cfg.fine_range = config.fuzzy_fine_range;
  fuzzy_cfg.diff = config.fuzzy_diff;
  fuzzy_cfg.n_bins = config.fuzzy_n_bins;
  auto fuzzy_quantizer = signal::make_fuzzy_quantizer(fuzzy_cfg);
  if (!fuzzy_quantizer) {
    throw std::runtime_error("Failed to create fuzzy quantizer: " + config.fuzzy_quantizer);
  }

  signal::SeedExtractorConfig extractor_cfg;
  extractor_cfg.backend = config.seed_type;
  extractor_cfg.k = config.seed_k;
  extractor_cfg.stride = config.seed_stride;
  extractor_cfg.window = config.minimizer_window;
  extractor_cfg.qbits = 4;
  auto extractor = signal::make_seed_extractor(extractor_cfg);
  if (!extractor) {
    throw std::runtime_error("Failed to create seed extractor");
  }

  IndexPipelineResult result;
  std::vector<std::size_t> path_lengths;

  BucketIndexConfig bi_config;
  bi_config.seed_k = config.seed_k;
  bi_config.seed_stride = config.seed_stride;
  bi_config.executor = config.executor;

  auto bi_result = bucketIndex(flat_graph, model, *fuzzy_quantizer, *extractor, bi_config);

  stage_elapsed =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - stage_start)
          .count();
  LOG_INFO("[2/2] indexed: " + std::to_string(bi_result.seed_store->size()) +
           " unique seeds [" + std::to_string(stage_elapsed) + "s]");

  path_lengths = std::move(bi_result.path_lengths);
  result.seed_store = std::move(bi_result.seed_store);
  result.linearization_coords = std::move(bi_result.linearization_coords);

  // Store path lengths in FlatGraph
  for (std::size_t i = 0; i < flat_graph.pathCount(); ++i) {
    flat_graph.setPathLength(static_cast<std::uint32_t>(i), path_lengths[i]);
  }

  // Package graph store
  result.graph_store = std::make_unique<FlatGraphStore>(std::move(flat_graph));

  /* 3. Compute 1D sort coordinates (for SortChainer) */

  if (config.compute_1d_sort) {
    auto sort_start = std::chrono::high_resolution_clock::now();
    result.node_1d_coords = compute_1d_sort(result.graph_store->flat(),
                                             result.linearization_coords,
                                             path_lengths, config.sort_1d_config);
    auto sort_elapsed =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - sort_start)
            .count();
    LOG_INFO("[3/3] 1D sort: " + std::to_string(result.node_1d_coords.size()) +
             " node positions [" + std::to_string(sort_elapsed) + "s]");
  }
  result.path_lengths = std::move(path_lengths);
  result.pore_k = model.k();
  result.model_name = model.name();
  result.fuzzy_quantizer = fuzzy_cfg.backend;

  return result;
}

}  // namespace piru::index
