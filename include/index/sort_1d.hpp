/**
 * sort_1d.hpp
 *
 * 1D graph linearization via path-guided stochastic gradient descent (PG-SGD).
 * Assigns each node a single 1D coordinate that approximates nucleotide
 * distances along paths. Used by SortChainer for fast approximate chaining.
 *
 * Algorithm based on odgi's path_linear_sgd (path_sgd.cpp) which implements
 * "Graph Drawing by Stochastic Gradient Descent" (Zheng et al. 2018,
 * https://arxiv.org/abs/1710.04626).
 *
 * Related:
 *  - sort_1d.cpp
 *  - linearizer.hpp (provides path coordinates used as target distances)
 *  - index_pipeline.hpp (calls compute_1d_sort)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "index/flat_graph.hpp"
#include "index/linearizer.hpp"

namespace piru::index {

struct Sort1DConfig {
  std::uint64_t iter_max{100};        // max SGD iterations
  std::uint64_t iter_with_max_lr{0};  // iteration with max learning rate
  double eps{0.01};                   // min learning rate factor
  double delta{0.0};                  // convergence threshold (0 = run all iters)
  double theta{0.99};                 // Zipfian skew parameter (matches odgi)
  std::uint64_t space{0};             // max sampling distance in path steps (0 = longest path)
  std::uint64_t space_max{100};       // threshold for zeta quantization
  double cooling_start{0.5};          // fraction of iterations before cooling phase
  std::uint64_t seed{9399220};        // RNG seed
};

/**
 * Compute 1D coordinates for all nodes via path-guided SGD.
 *
 * @param graph       The alignment graph (for node count, sequence lengths)
 * @param coords      Per-node linearization coordinates (for path distances)
 * @param path_lengths Length of each path
 * @param config      SGD parameters
 * @return            node_1d_coord[node_id] = 1D position (double)
 */
std::vector<float> compute_1d_sort(const FlatGraph& graph,
                                   const std::vector<std::vector<LinearCoordinate>>& coords,
                                   const std::vector<std::size_t>& path_lengths,
                                   const Sort1DConfig& config = {});

/**
 * Import 1D coordinates from odgi layout TSV.
 * odgi format: idx\tX\tY (two rows per node: start and end positions, Y=0).
 * Maps odgi's original node ordering to our directional graph (fwd/rev pairs).
 *
 * @param path        TSV file path
 * @param num_nodes   Total node count in directional graph (2 * original nodes)
 * @return            node_1d_coord[node_id] for all nodes (fwd and rev)
 */
std::vector<float> import_1d_coords_odgi(const std::string& path, std::size_t num_nodes);

/**
 * Dump 1D coordinates to TSV for visualization/validation.
 * Format: node_id\tstart_pos\tend_pos\tseq_len
 */
void dump_1d_coords_tsv(const std::string& path, const std::vector<float>& coords,
                        const FlatGraph& graph);

}  // namespace piru::index
