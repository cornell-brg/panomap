/**
 * graph_chainer2.hpp
 *
 * Graph-space DP chainer v2: per-path sorted iteration with shared DP table.
 * Equivalent to PathChainer but chains NodeAnchors (no alias chains) and
 * allows haplotype hopping via the shared DP table.
 *
 * Algorithm:
 * 1. Sort NodeAnchors by read_pos (DP ordering)
 * 2. Build per-path sorted lists (by ref_coord) with back-pointers to dp_idx
 * 3. For each anchor i (read_pos order), for each path p it's on:
 *    - Find anchor i's position in path p's sorted list
 *    - Scan backwards in path p's list for predecessors within band
 *    - Update shared dp[i] with best score across all paths
 * 4. Extract chains with truncate-on-used backtracking
 *
 * Related:
 *  - graph_chainer.hpp  (v1: binary search approach)
 *  - path_chainer.hpp   (per-path chainer this design is modeled after)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <vector>

#include "index/linearizer.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

class GraphChainer2 : public Chainer {
public:
  GraphChainer2(const std::vector<std::vector<index::LinearCoordinate>>& coords,
                const std::vector<std::size_t>& path_lengths,
                std::size_t max_dist = 500, std::size_t max_diag_dev = 500,
                std::size_t min_chain_score = 12, std::size_t max_chains = 10,
                std::size_t max_skip = 25);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "graph-chain2"; }

private:
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;
  index::TransposedCoords transposed_;

  std::size_t max_dist_;
  std::size_t max_diag_dev_;
  std::size_t min_chain_score_;
  std::size_t max_chains_;
  std::size_t max_skip_;

  struct DPEntry {
    double score{0.0};
    int pred{-1};
    std::size_t path_id{0};
    std::int64_t ref_coord{0};
  };

  double anchorScore(const NodeAnchor& anchor) const;
};

}  // namespace piru::mapping
