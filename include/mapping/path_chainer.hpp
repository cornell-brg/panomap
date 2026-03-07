// SPDX-License-Identifier: MIT
// Path-space DP chainer: expands to per-path anchors, chains per path.

#pragma once

#include <cstddef>
#include <vector>

#include "cli/parse.hpp"
#include "index/linearizer.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

// Configuration for path-space chaining algorithm.
// Defaults tuned for noisy nanopore signals (DEV027).
struct PathChainerConfig {
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
  static PathChainerConfig from_parsed(const cli::Parsed& parsed);
};

// Grouped PathAnchors by path_id. path_id is implicit (the vector index).
using PathAnchorGroups = std::vector<std::vector<PathAnchor>>;

// Path-space colinear chainer (Method 1).
//
// Internally: expand NodeAnchors to PathAnchors grouped by path_id,
// optionally merge adjacent anchors, then DP chain per path.
//
// Algorithm:
// 1. Expand NodeAnchors -> PathAnchorGroups (grouped by path_id)
// 2. Merge adjacent/overlapping PathAnchors per path (optional)
// 3. Per path: sort by (ref_coord, query_pos), DP chain
// 4. Collect best chains across all paths, rank by score
class PathChainer : public Chainer {
public:
  // coords[node_id] = linearization coordinates for that node (non-owning)
  // path_lengths[path_id] = length of that path for bounds checking (non-owning)
  PathChainer(PathChainerConfig config, const std::vector<std::vector<index::LinearCoordinate>>& coords,
            const std::vector<std::size_t>& path_lengths);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "path-chain"; }

private:
  PathChainerConfig config_;
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;

  // Expand NodeAnchors to PathAnchors grouped by path_id
  PathAnchorGroups expand(const std::vector<NodeAnchor>& hits) const;

  // DP internals (operate on per-path PathAnchors)
  // hits passed for back-reference via src_idx during chain extraction
  ChainResult chain_grouped(const PathAnchorGroups& groups, const std::vector<NodeAnchor>& hits) const;
  std::vector<Chain> chain_one_path(const std::vector<PathAnchor>& anchors, std::size_t path_id,
                                    const std::vector<NodeAnchor>& hits) const;
  bool can_chain(const PathAnchor& j, const PathAnchor& i) const;
  double gap_cost(const PathAnchor& j, const PathAnchor& i) const;
  double anchor_score(const PathAnchor& anchor) const;
  std::vector<std::size_t> backtrack_chain(const std::vector<int>& pred,
                                           std::size_t best_idx) const;
};

}  // namespace piru::mapping
