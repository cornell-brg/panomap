// SPDX-License-Identifier: MIT
// Unified mapping result types for the mapping pipeline.

#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "alignment/segment_aligner.hpp"  // GraphPosition
#include "mapping/seed_clusterer.hpp"     // SeedAnchor

namespace piru::mapping {

/// A single mapping (one chain, optionally aligned).
///
/// After DP chaining: anchors and chain_score are populated.
/// After alignment (if enabled): alignment fields are populated.
struct Mapping {
  // Always present (from DP chaining)
  std::vector<SeedAnchor> anchors;
  double chain_score{0.0};

  // Optional (present if --align enabled)
  std::optional<float> alignment_cost;
  std::optional<float> normalized_alignment_cost;
  std::optional<std::vector<alignment::GraphPosition>> alignment_path;
  std::size_t segments_aligned{0};

  /// Check if alignment was performed.
  bool hasAlignment() const { return alignment_cost.has_value(); }
};

/// All mappings for a single read.
struct ReadMapResult {
  std::vector<Mapping> mappings;          // primary (index 0) + secondaries
  std::size_t expanded_anchor_count{0};   // anchors before chaining

  /// Check if read mapped (has at least one mapping).
  bool mapped() const { return !mappings.empty(); }

  /// Get primary mapping (first one, if any).
  const Mapping* primary() const {
    return mappings.empty() ? nullptr : &mappings[0];
  }
};

}  // namespace piru::mapping
