/**
 * sort_1d.cpp
 *
 * Path-guided SGD for 1D graph linearization. Closely follows odgi's
 * path_linear_sgd (related/odgi/src/algorithms/path_sgd.cpp).
 *
 * The algorithm iteratively adjusts 1D node positions so that the distance
 * between node pairs in the layout matches their nucleotide distance along
 * shared paths. Uses Zipfian sampling to prioritize nearby node pairs,
 * with 50% uniform sampling for long-range corrections.
 *
 * Related:
 *  - sort_1d.hpp
 *  - third_party/dirtyzipf/dirty_zipfian_int_distribution.h
 *  - "Graph Drawing by Stochastic Gradient Descent" (Zheng et al. 2018)
 *
 * SPDX-License-Identifier: MIT
 */

#include "index/sort_1d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

#include "dirty_zipfian_int_distribution.h"
#include "util/logging.hpp"

namespace piru::index {

namespace {

// Build a path step index for SGD sampling.
// For each path, store the ordered list of (node_id, cumulative_ref_pos).
struct PathStep {
  std::uint32_t node_id;
  std::int64_t ref_pos;  // cumulative position along path
};
using PathStepIndex = std::vector<std::vector<PathStep>>;

PathStepIndex build_path_step_index(
    const AlnGraph& /*graph*/,
    const std::vector<std::vector<LinearCoordinate>>& coords) {
  // Determine number of paths
  std::size_t num_paths = 0;
  for (const auto& node_coords : coords) {
    for (const auto& lc : node_coords) {
      num_paths = std::max(num_paths, lc.path_id + 1);
    }
  }

  // Collect steps per path, using only forward nodes (even node IDs)
  // and forward paths (even path IDs). simpleExpand creates both forward
  // and reverse versions; we only need forward to match odgi's behavior.
  PathStepIndex index(num_paths);
  for (std::size_t nid = 0; nid < coords.size(); nid += 2) {
    std::uint32_t fwd_idx = static_cast<std::uint32_t>(nid / 2);
    for (const auto& lc : coords[nid]) {
      if (lc.path_id % 2 != 0) continue;  // skip reverse-strand paths
      index[lc.path_id].push_back({fwd_idx, lc.ref_coord});
    }
  }

  // Sort each path's steps by ref_pos
  for (auto& steps : index) {
    std::sort(steps.begin(), steps.end(),
              [](const PathStep& a, const PathStep& b) {
                return a.ref_pos < b.ref_pos;
              });
  }

  return index;
}

// Compute learning rate schedule (matches odgi's path_linear_sgd_schedule).
std::vector<double> compute_schedule(double w_min, double w_max,
                                     std::uint64_t iter_max,
                                     std::uint64_t iter_with_max_lr,
                                     double eps) {
  double eta_max = 1.0 / w_min;
  double eta_min = eps / w_max;
  double lambda = std::log(eta_max / eta_min) / static_cast<double>(iter_max - 1);

  std::vector<double> etas;
  etas.reserve(iter_max + 1);
  for (std::uint64_t t = 0; t <= iter_max; ++t) {
    double exp_arg = -lambda * std::abs(static_cast<int64_t>(t) -
                                        static_cast<int64_t>(iter_with_max_lr));
    etas.push_back(eta_max * std::exp(exp_arg));
  }
  return etas;
}

}  // namespace

std::vector<float> compute_1d_sort(
    const AlnGraph& graph,
    const std::vector<std::vector<LinearCoordinate>>& coords,
    const std::vector<std::size_t>& /*path_lengths*/,
    const Sort1DConfig& config) {
  const std::size_t num_nodes = graph.nodeCount();

  /* 1. Initialize positions: cumulative sequence length (matches odgi).
   * Only forward nodes participate in SGD. Reverse nodes derive positions
   * from their forward partner at the end.
   * Node ID scheme from simpleExpand: forward = orig*2, reverse = orig*2+1. */
  std::size_t num_fwd = num_nodes / 2;
  std::vector<double> X(num_fwd);
  std::vector<std::size_t> fwd_lengths(num_fwd);
  std::uint64_t cumulative = 0;
  for (std::size_t i = 0; i < num_fwd; ++i) {
    X[i] = static_cast<double>(cumulative);
    fwd_lengths[i] = graph.node(i * 2).sequence.size();
    cumulative += fwd_lengths[i];
  }

  /* 2. Build path step index for sampling. */
  auto path_index = build_path_step_index(graph, coords);

  // Filter to paths with more than one step
  std::vector<std::size_t> valid_paths;
  for (std::size_t p = 0; p < path_index.size(); ++p) {
    if (path_index[p].size() > 1) {
      valid_paths.push_back(p);
    }
  }
  if (valid_paths.empty()) {
    LOG_WARN("No valid paths for 1D SGD (all paths have <= 1 step)");
    return std::vector<float>(X.begin(), X.end());
  }

  // Total steps across all valid paths (for uniform sampling)
  std::uint64_t total_steps = 0;
  std::size_t max_path_steps = 0;
  for (auto p : valid_paths) {
    total_steps += path_index[p].size();
    max_path_steps = std::max(max_path_steps, path_index[p].size());
  }

  /* 3. Determine space (max Zipfian jump distance in path steps).
   * Default (space=0): longest path step count, matching odgi. */
  std::uint64_t space = config.space;
  if (space == 0) space = max_path_steps;
  std::uint64_t space_max = config.space_max;

  /* 4. Precompute zeta values for Zipfian distributions.
   * Matches odgi: precompute for all jump distances up to space,
   * with quantization beyond space_max. */
  std::uint64_t num_zetas = (space <= space_max)
      ? space
      : space_max + (space - space_max) / space_max + 1;
  std::vector<double> zetas(num_zetas + 1, 0.0);
  {
    double zeta_tmp = 0.0;
    for (std::uint64_t i = 1; i <= space; ++i) {
      zeta_tmp += dirtyzipf::fast_precise_pow(1.0 / static_cast<double>(i), config.theta);
      if (i <= space_max) {
        zetas[i] = zeta_tmp;
      }
      if (i >= space_max && (i - space_max) % space_max == 0) {
        zetas[space_max + 1 + (i - space_max) / space_max] = zeta_tmp;
      }
    }
  }

  /* 5. Compute learning rate schedule (matches odgi defaults). */
  double eta_max_val = static_cast<double>(max_path_steps) * static_cast<double>(max_path_steps);
  double w_min = 1.0 / eta_max_val;
  double w_max = 1.0;
  auto etas = compute_schedule(w_min, w_max, config.iter_max,
                                config.iter_with_max_lr, config.eps);

  /* 6. SGD loop */
  std::uint64_t min_term_updates = total_steps;
  std::uint64_t first_cooling_iter = static_cast<std::uint64_t>(
      std::floor(config.cooling_start * static_cast<double>(config.iter_max)));

  std::mt19937_64 rng(config.seed);

  // Build cumulative step count for weighted path sampling
  std::vector<std::uint64_t> cum_steps;
  cum_steps.reserve(valid_paths.size());
  std::uint64_t running = 0;
  for (auto p : valid_paths) {
    running += path_index[p].size();
    cum_steps.push_back(running);
  }
  std::uniform_int_distribution<std::uint64_t> global_step_dist(0, running - 1);
  std::uniform_int_distribution<int> flip(0, 1);  // coin flip for 50/50 sampling

  double delta_max = 0.0;

  for (std::uint64_t iteration = 0; iteration < config.iter_max; ++iteration) {
    double eta = etas[iteration];
    delta_max = 0.0;

    // Adjust theta for cooling phase
    bool cooling = iteration > first_cooling_iter;
    double adj_theta = cooling ? 0.001 : config.theta;

    for (std::uint64_t update = 0; update < min_term_updates; ++update) {
      /* Sample first step uniformly from all path steps */
      std::uint64_t global_idx = global_step_dist(rng);
      auto it = std::lower_bound(cum_steps.begin(), cum_steps.end(), global_idx + 1);
      std::size_t path_rank = static_cast<std::size_t>(it - cum_steps.begin());
      std::size_t p = valid_paths[path_rank];
      const auto& steps = path_index[p];
      std::size_t step_count = steps.size();
      std::size_t offset_in_path = (path_rank > 0) ? cum_steps[path_rank - 1] : 0;
      std::size_t s_rank = static_cast<std::size_t>(global_idx - offset_in_path);

      /* Sample second step: 50% Zipfian (nearby) + 50% uniform (long-range).
       * Matches odgi: if (cooling || flip(gen)) -> Zipfian, else -> uniform.
       * During cooling, always Zipfian with theta=0.001 (near-uniform). */
      std::size_t s_rank_b;
      if (cooling || flip(rng)) {
        // Zipfian sampling: bias toward nearby steps
        std::uint64_t jump_space;
        bool go_backward = (s_rank > 0 && flip(rng)) || s_rank == step_count - 1;

        if (go_backward) {
          jump_space = std::min(space, static_cast<std::uint64_t>(s_rank));
        } else {
          jump_space = std::min(space, static_cast<std::uint64_t>(step_count - s_rank - 1));
        }
        if (jump_space == 0) continue;

        // Look up precomputed zeta (with quantization for large jump_space)
        std::uint64_t zeta_idx = jump_space;
        if (jump_space > space_max) {
          zeta_idx = space_max + 1 + (jump_space - space_max) / space_max;
        }
        if (zeta_idx >= zetas.size()) zeta_idx = zetas.size() - 1;

        dirtyzipf::dirty_zipfian_int_distribution<std::uint64_t>::param_type z_p(
            1, jump_space, adj_theta, zetas[zeta_idx]);
        dirtyzipf::dirty_zipfian_int_distribution<std::uint64_t> z(z_p);
        std::uint64_t z_i = z(rng);

        if (go_backward) {
          s_rank_b = s_rank - static_cast<std::size_t>(z_i);
        } else {
          s_rank_b = s_rank + static_cast<std::size_t>(z_i);
        }
      } else {
        // Uniform sampling: random step anywhere on this path (long-range correction)
        std::uniform_int_distribution<std::size_t> uni_dist(0, step_count - 1);
        s_rank_b = uni_dist(rng);
      }

      if (s_rank_b >= step_count || s_rank_b == s_rank) continue;

      // Get nodes and target distance
      std::uint32_t node_i = steps[s_rank].node_id;
      std::uint32_t node_j = steps[s_rank_b].node_id;
      if (node_i == node_j) continue;

      double target_dist = std::abs(
          static_cast<double>(steps[s_rank].ref_pos) -
          static_cast<double>(steps[s_rank_b].ref_pos));
      if (target_dist == 0) continue;

      // Current distance in layout
      double dx = X[node_i] - X[node_j];
      if (dx == 0) dx = 1e-9;
      double mag = std::abs(dx);

      // Weight inversely proportional to distance (matches odgi)
      double w_ij = 1.0 / target_dist;
      double mu = eta * w_ij;
      if (mu > 1.0) mu = 1.0;

      // SGD update
      double delta = mu * (mag - target_dist) / 2.0;
      double r = delta / mag;
      double r_x = r * dx;

      X[node_i] -= r_x;
      X[node_j] += r_x;

      double delta_abs = std::abs(delta);
      if (delta_abs > delta_max) delta_max = delta_abs;
    }

    // Convergence check
    if (config.delta > 0.0 && delta_max <= config.delta) {
      LOG_INFO("1D SGD converged at iteration " + std::to_string(iteration) +
               " (delta_max=" + std::to_string(delta_max) + ")");
      break;
    }
  }

  LOG_INFO("1D SGD complete: " + std::to_string(num_fwd) + " forward nodes, " +
           std::to_string(config.iter_max) + " iterations");

  // Expand to full node count and convert double -> float32
  std::vector<float> X_full(num_nodes);
  for (std::size_t i = 0; i < num_fwd; ++i) {
    X_full[i * 2] = static_cast<float>(X[i]);                                   // forward: start
    X_full[i * 2 + 1] = static_cast<float>(X[i] + static_cast<double>(fwd_lengths[i]));  // reverse: end
  }

  return X_full;
}

std::vector<float> import_1d_coords_odgi(const std::string& path,
                                          std::size_t num_nodes) {
  // odgi TSV: idx\tX\tY, two rows per node (start, end). Y=0 for 1D layout.
  // odgi has N original nodes, we have 2N directional nodes.
  std::size_t num_orig = num_nodes / 2;
  std::vector<float> result(num_nodes, 0.0f);

  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open 1D coords file: " + path);
  }

  std::string header;
  std::getline(ifs, header);  // skip "idx\tX\tY"

  std::size_t orig_idx = 0;
  std::string line;
  while (std::getline(ifs, line) && orig_idx < num_orig) {
    // Read start line
    double x_start = 0.0;
    {
      auto tab1 = line.find('\t');
      auto tab2 = line.find('\t', tab1 + 1);
      x_start = std::stod(line.substr(tab1 + 1, tab2 - tab1 - 1));
    }

    // Read end line
    double x_end = x_start;
    if (std::getline(ifs, line)) {
      auto tab1 = line.find('\t');
      auto tab2 = line.find('\t', tab1 + 1);
      x_end = std::stod(line.substr(tab1 + 1, tab2 - tab1 - 1));
    }

    // Forward node: start position. Reverse node: end position.
    result[orig_idx * 2] = static_cast<float>(x_start);
    result[orig_idx * 2 + 1] = static_cast<float>(x_end);
    ++orig_idx;
  }

  LOG_INFO("Imported 1D coords: " + std::to_string(orig_idx) +
           " original nodes -> " + std::to_string(num_nodes) + " directional nodes");

  return result;
}

void dump_1d_coords_tsv(const std::string& path,
                         const std::vector<float>& coords,
                         const AlnGraph& graph) {
  std::ofstream ofs(path);
  ofs << "node_id\tstart_pos\tend_pos\tseq_len\n";
  for (std::size_t i = 0; i < coords.size(); ++i) {
    std::size_t seq_len = graph.node(i).sequence.size();
    ofs << i << "\t" << coords[i] << "\t" << (coords[i] + seq_len) << "\t" << seq_len << "\n";
  }
}

}  // namespace piru::index
