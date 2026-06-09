/**
 * path_chainer.hpp
 *
 * Path-space chainer, per-path DP, component-aware interval dedup.
 *
 * Related:
 *  - path_chainer.cpp
 *  - chainer.hpp
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "cli/parse.hpp"
#include "core/index/linearizer.hpp"
#include "core/mapping/chainer.hpp"

namespace panomap::mapping {

// Configuration for path-space chaining algorithm.
// Scoring aligned with RawHash2/minimap2 (DEV061).
struct PathChainerConfig {
  std::size_t max_dist_ref{2500};      // Max ref distance between anchors
  std::size_t max_dist_query{2500};    // Max query distance between anchors
  std::size_t bw{500};                 // Max diagonal deviation |dr - dq| (bandwidth)
  std::size_t min_chain_score{15};     // Min score to report a chain
  std::size_t min_chain_anchors{2};    // Min anchors per chain
  std::size_t max_chains{10};          // Hard safety cap (internal, not a CLI knob)
  float secondary_ratio{0.3f};        // Stop extracting when score < ratio * primary
  std::size_t max_survivor_chains{0};  // Max chains for survivor marking (0 = unlimited)
  std::size_t max_skip{25};            // Stop after N consecutive failed chain attempts
  std::size_t max_iterations{0};       // Max predecessors to check per anchor (0 = unlimited)

  float chn_pen_gap{0.8f};   // Gap penalty factor (applied to diagonal deviation)
  float chn_pen_skip{0.0f};  // Skip penalty factor (applied to gap distance)

  bool merge_chains{false};  // Merge overlapping chains on same path
  bool merge_anchors{true};  // Merge overlapping anchors before chaining
  std::size_t pore_k{0};     // Pore model k. When >0, scoring span = pore_k + anchor_len - 1.

  // CLI integration: options and parsing for --chain-* flags.
  static std::vector<cli::Option> cli_options();
  static PathChainerConfig from_parsed(const cli::Parsed& parsed);
};

// Path-space colinear chainer (Method 1).
//
// Internally: expand NodeAnchors to PathAnchors grouped by path_id,
// optionally merge adjacent anchors, then DP chain per path.
//
// Algorithm:
// 1. Expand NodeAnchors -> PathAnchorGroups (grouped by path_id)
// 2. Merge adjacent/overlapping PathAnchors per path (optional)
// 3. Per path: sort by (ref_coord, query_pos), DP chain (SoA layout)
// 4. Collect best chains across all paths, rank by score
class PathChainer : public Chainer {
public:
  // coords[node_id] = linearization coordinates for that node (non-owning)
  // path_lengths[path_id] = length of that path for bounds checking (non-owning)
  // node_1d_coords, component_ids: optional, enables interval-overlap dedup
  PathChainer(PathChainerConfig config,
              const std::vector<std::vector<index::LinearCoordinate>>& coords,
              const std::vector<std::size_t>& path_lengths,
              const std::vector<float>* node_1d_coords = nullptr,
              const std::vector<std::uint32_t>* component_ids = nullptr);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "path-chain"; }

private:
  PathChainerConfig config_;
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;
  const std::vector<float>* node_1d_coords_;
  const std::vector<std::uint32_t>* component_ids_;
};

}  // namespace panomap::mapping
