// SPDX-License-Identifier: MIT
// Formats mapping results to PAF/GAF/GAM output.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "io/results/result.hpp"
#include "mapping/map_result.hpp"

namespace piru::mapping {

/// Configuration for result formatting.
struct ResultFormatterConfig {
  bool primary_only{false};      // Only output primary (best) mapping
  std::size_t max_secondary{9};  // Max secondary mappings (default: 9, so up to 10 total)
};

/// Formats ReadMapResult to AlignmentResult for PAF/GAF output.
///
/// Pure data transformation - no computation, just field mapping.
/// Alignment is already done in process_read() if enabled.
class ResultFormatter {
 public:
  explicit ResultFormatter(const index::AlnGraph& graph,
                           ResultFormatterConfig config = {});

  /// Format mappings to AlignmentResult vector.
  /// Returns one result per mapping (up to max_secondary + 1).
  std::vector<io::AlignmentResult> format(
      const ReadMapResult& map_result,
      const std::string& read_id,
      std::size_t read_length) const;

 private:
  /// Format a single Mapping to AlignmentResult.
  io::AlignmentResult formatMapping(
      const Mapping& mapping,
      const std::string& read_id,
      std::size_t read_length,
      bool is_primary) const;

  /// Build GAF-style path string from anchors (">node1>node2>node3").
  std::string buildPathString(const std::vector<SeedAnchor>& anchors) const;

  /// Get path name from path_id.
  const std::string& getPathName(std::size_t path_id) const;

  /// Compute path length (sum of node lengths along path).
  std::size_t computePathLength(std::size_t path_id) const;

  /// Estimate MAPQ from chain scores.
  int estimateMapQ(double primary_score, double secondary_score) const;

  const index::AlnGraph& graph_;
  ResultFormatterConfig config_;

  // Cached path lengths (computed on first access).
  mutable std::vector<std::size_t> path_lengths_;
  mutable bool path_lengths_computed_{false};
};

}  // namespace piru::mapping
