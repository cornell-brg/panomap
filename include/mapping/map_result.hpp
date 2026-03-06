// SPDX-License-Identifier: MIT
// Unified mapping result types for the mapping pipeline.

#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "mapping/chainer.hpp"  // ChainedAnchor

namespace piru::mapping {

/// A single mapping (one chain, optionally aligned).
///
/// After DP chaining: anchors and chain_score are populated.
/// After alignment (if enabled): alignment fields are populated.
struct Mapping {
  // Always present (from DP chaining)
  std::vector<ChainedAnchor> anchors;
  double chain_score{0.0};

  // Optional alignment scores (reserved for future use)
  std::optional<float> alignment_cost;
  std::optional<float> normalized_alignment_cost;
};

/// All mappings for a single read.
struct ReadMapResult {
  std::vector<Mapping> mappings;         // primary (index 0) + secondaries
  std::size_t expanded_anchor_count{0};  // anchors before chaining

  // ROI classification (populated when --roi is active)
  double roi_overlap{-1.0};  // -1 = not computed
  bool roi_keep{false};      // final keep/reject decision

  /// Check if read mapped (has at least one mapping).
  bool mapped() const { return !mappings.empty(); }

  /// Get primary mapping (first one, if any).
  const Mapping* primary() const { return mappings.empty() ? nullptr : &mappings[0]; }
};

}  // namespace piru::mapping
