/**
 * graph_chainer.hpp
 *
 * Graph-space chainer with haplotype hopping.
 *
 * Related:
 *  - graph_chainer.cpp
 *  - chainer.hpp
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <vector>

#include "index/linearizer.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

// Graph-space colinear chainer (Method 2).
//
// Chains NodeAnchors directly in graph space. For each anchor pair (j, i),
// finds the best (path, ref_coord) pair by querying the coord index.
// Eliminates alias chains and enables haplotype hopping.
//
// Algorithm:
// 1. Sort NodeAnchors by read_pos
// 2. DP: for each pair (j, i), iterate j's path coords, probe i on same
//    path, compute gap cost. Band on query and ref distance.
// 3. Track best score, predecessor, and chosen path per anchor
// 4. Extract chains, rank by score
class GraphChainer : public Chainer {
public:
  // coords[node_id] = linearization coordinates for that node (non-owning)
  // path_lengths[path_id] = length of that path for bounds checking (non-owning)
  GraphChainer(const std::vector<std::vector<index::LinearCoordinate>>& coords,
               const std::vector<std::size_t>& path_lengths, std::size_t max_dist = 500,
               std::size_t max_diag_dev = 500, std::size_t min_chain_score = 12,
               std::size_t max_chains = 10, std::size_t max_skip = 25);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "graph-chain"; }

private:
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;
  index::TransposedCoords transposed_;

  std::size_t max_dist_;
  std::size_t max_diag_dev_;
  std::size_t min_chain_score_;
  std::size_t max_chains_;
  std::size_t max_skip_;

  /* DP state per anchor */
  struct DPEntry {
    double score{0.0};
    int pred{-1};               // predecessor index (-1 = chain start)
    std::size_t path_id{0};     // chosen path for this anchor
    std::int64_t ref_coord{0};  // chosen ref_coord on that path
  };

  /* Find best transition from anchor j to anchor i across all shared paths.
   * Returns (gap_cost, path_id, ref_coord_j, ref_coord_i) or negative if no
   * valid transition exists. */
  struct Transition {
    double cost;
    std::size_t path_id;
    std::int64_t ref_j;
    std::int64_t ref_i;
    bool valid{false};
  };
  Transition bestTransition(const NodeAnchor& j, const NodeAnchor& i,
                            std::size_t prev_path_id) const;

  double anchorScore(const NodeAnchor& anchor) const;
  std::vector<std::size_t> backtrack(const std::vector<DPEntry>& dp, std::size_t best_idx) const;
};

}  // namespace piru::mapping
