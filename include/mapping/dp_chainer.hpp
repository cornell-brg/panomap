// SPDX-License-Identifier: MIT
// DP-based colinear chainer for path-aware seed selection.

#pragma once

#include <cstddef>
#include <vector>

#include "cli/parse.hpp"
#include "index/linearizer.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

// Configuration for DP-based chaining algorithm.
// Defaults tuned for noisy nanopore signals (DEV027).
struct DPChainerConfig {
  std::size_t max_dist{500};        // Max query/ref distance for chaining (banding)
  std::size_t max_diag_dev{500};    // Max diagonal deviation |dr - dq|
  std::size_t min_chain_score{12};  // Min score to report a chain
  std::size_t max_chains{10};       // Max number of chains to extract (multi-mapping)
  std::size_t max_skip{25};         // Stop after N consecutive failed chain attempts

  double anchor_weight{1.0};            // Weight per anchor length
  double gap_penalty_factor{0.02};      // Penalty per unit gap distance
  double diag_penalty_factor{0.05};     // Penalty per unit diagonal deviation
  double overlap_penalty_factor{0.90};  // Penalty per unit overlap

  bool merge_chains{false};  // Merge overlapping chains on same path
  bool merge_anchors{true};  // Merge overlapping anchors before chaining

  // CLI integration: options and parsing for --chain-* flags.
  static std::vector<cli::Option> cli_options();
  static DPChainerConfig from_parsed(const cli::Parsed& parsed);
};

// DP-based colinear chainer.
//
// Internally: expand NodeAnchors to PathAnchors (linear space), optionally
// merge adjacent anchors, then DP chain.
//
// Algorithm:
// 1. Expand NodeAnchors -> PathAnchors using linearization coordinates
// 2. Merge adjacent/overlapping PathAnchors (optional)
// 3. Sort by (path_id, ref_coord, query_pos)
// 4. DP: dp[i] = max_j(dp[j] + score(i) - gap_cost(j, i))
// 5. Backtrack to extract chain(s)
class DPChainer : public Chainer {
public:
  // coords[node_id] = linearization coordinates for that node (non-owning)
  // path_lengths[path_id] = length of that path for bounds checking (non-owning)
  DPChainer(DPChainerConfig config, const std::vector<std::vector<index::LinearCoordinate>>& coords,
            const std::vector<std::size_t>& path_lengths);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "dp-chain"; }

private:
  DPChainerConfig config_;
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;

  // Expand NodeAnchors to PathAnchors using linearization coordinates
  std::vector<PathAnchor> expand(const std::vector<NodeAnchor>& hits) const;

  // DP internals (operate on PathAnchors after expansion)
  ChainResult chain_path_anchors(const std::vector<PathAnchor>& anchors) const;
  bool can_chain(const PathAnchor& j, const PathAnchor& i) const;
  double gap_cost(const PathAnchor& j, const PathAnchor& i) const;
  double anchor_score(const PathAnchor& anchor) const;
  std::vector<std::size_t> backtrack_chain(const std::vector<int>& pred,
                                           std::size_t best_idx) const;
};

}  // namespace piru::mapping
