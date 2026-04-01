/**
 * pan_chainer.hpp
 *
 * PanChainer: 1D-banded cross-path colinear chaining. Sorts anchors by
 * canonical 1D coord (PG-SGD), bands in 1D for cheap candidate selection,
 * then evaluates colinearity using exact per-path ref coords.
 *
 * Combines SortChainer's O(anchors) scaling with PathChainer's accuracy.
 * Single DP value per anchor -- path choice is local to each transition.
 *
 * Related:
 *  - pan_chainer.cpp
 *  - sort_chainer.hpp (1D-only chaining)
 *  - path_chainer.hpp (per-path chaining)
 *  - knowledge/scratch/note-260325-panchainer-dp-formulation.md
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "cli/parse.hpp"
#include "index/linearizer.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

struct PanChainerConfig {
  std::size_t band_1d{5000};         // 1D band width for candidate selection
  std::size_t max_dist_ref{2500};    // max ref distance on shared path
  std::size_t max_dist_query{2500};  // max query distance between anchors
  std::size_t bw{500};               // max diagonal deviation |dr - dq|
  std::size_t min_chain_score{15};
  std::size_t min_chain_anchors{2};
  std::size_t max_chains{10};
  std::size_t max_skip{25};
  std::size_t max_iterations{0};

  float chn_pen_gap{0.8f};  // same as PathChainer (exact path coords)
  float chn_pen_skip{0.0f};
  float chn_pen_switch{50.0f};  // penalty for switching haplotype path between transitions

  std::size_t pore_k{0};

  static PanChainerConfig from_parsed(const cli::Parsed& parsed);
};

/**
 * PanChainer: 1D-banded cross-path colinear chaining.
 *
 * Sort by canonical 1D coord, band in 1D, score using exact per-path
 * ref coords. Allows haplotype hopping -- consecutive chain links may
 * use different shared paths.
 */
class PanChainer : public Chainer {
public:
  PanChainer(PanChainerConfig config, const std::vector<float>& node_1d_coords,
             const std::vector<std::vector<index::LinearCoordinate>>& linearization_coords,
             const std::vector<std::size_t>& path_lengths);

  ChainResult chain(const std::vector<NodeAnchor>& hits) const override;
  std::string name() const override { return "pan-chain"; }

private:
  PanChainerConfig config_;
  const std::vector<float>& node_1d_coords_;
  const std::vector<std::vector<index::LinearCoordinate>>& coords_;
  const std::vector<std::size_t>& path_lengths_;
};

}  // namespace piru::mapping
