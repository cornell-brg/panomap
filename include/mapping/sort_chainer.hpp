/**
 * sort_chainer.hpp
 *
 * 1D sort-based chainer: chains NodeAnchors using pre-computed 1D
 * coordinates. Single DP pass, no per-path expansion. O(anchors) time.
 * Implicit haplotype hopping via the 1D graph linearization.
 *
 * "Fast mode" for adaptive sampling: maps/unmaps decision only,
 * no path_id assignment.
 *
 * Related:
 *  - sort_chainer.cpp
 *  - path_chainer.hpp (reference DP implementation)
 *  - sort_1d.hpp (1D coordinate computation)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <vector>

#include "cli/parse.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

struct SortChainerConfig {
  std::size_t max_dist_ref{5000};     // wider than PathChainer (2500) for 1D distortion
  std::size_t max_dist_query{5000};
  std::size_t bw{1000};              // wider bandwidth for 1D distortion
  std::size_t min_chain_score{15};
  std::size_t min_chain_anchors{2};
  std::size_t max_chains{10};
  std::size_t max_skip{25};

  float chn_pen_gap{0.4f};           // lower than PathChainer (0.8) for 1D distortion
  float chn_pen_skip{0.0f};
  float dd_tolerance_frac{0.0f};     // dead zone: dd below tolerance_frac * dg is penalty-free
  float chn_pen_ratio{0.5f};         // penalty for ratio inconsistency between consecutive transitions

  std::size_t pore_k{0};

  static SortChainerConfig from_parsed(const cli::Parsed& parsed);
};

/**
 * 1D sort-based chainer.
 *
 * Uses pre-computed 1D node coordinates to convert NodeAnchors to
 * (1d_ref_coord, read_pos) pairs. Single colinear DP pass.
 */
class SortChainer : public Chainer {
 public:
  SortChainer(SortChainerConfig config, const std::vector<float>& node_1d_coords,
              std::vector<std::uint32_t> node_bp_lens);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "sort-chain"; }

 private:
  SortChainerConfig config_;
  const std::vector<float>& node_1d_coords_;
  std::vector<std::uint32_t> node_bp_lens_;
};

}  // namespace piru::mapping
