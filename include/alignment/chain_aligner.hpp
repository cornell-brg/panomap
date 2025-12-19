// SPDX-License-Identifier: MIT
// ChainAligner: orchestrates anchor-to-anchor alignment for a full chain.

#pragma once

#include <memory>
#include <vector>

#include "alignment/segment_aligner.hpp"
#include "alignment/signal_utils.hpp"

namespace piru::alignment {

/// Result of aligning an entire chain.
struct ChainAlignmentResult {
  float total_cost;                    // Sum of segment costs
  std::vector<GraphPosition> path;     // Full path through graph
  size_t segments_aligned;             // Number of segments successfully aligned
  size_t segments_reached_target;      // Number of segments that reached target

  /// Check if alignment was successful (at least some segments aligned).
  bool valid() const { return segments_aligned > 0 && total_cost < std::numeric_limits<float>::infinity(); }

  /// Compute normalized cost (per query position).
  float normalizedCost(size_t query_length) const {
    return query_length > 0 ? total_cost / static_cast<float>(query_length) : 0.0f;
  }
};

/// Backend selection for segment alignment.
enum class AlignerBackend {
  kPathGuided,  // Linear extraction + standard DTW (fast, same-path only)
  kRadius,      // BFS extraction + graph DTW (handles bubbles)
  kAuto         // PathGuided for same-node, Radius for cross-node
};

/// Configuration for ChainAligner.
struct ChainAlignerConfig {
  AlignerBackend backend = AlignerBackend::kPathGuided;
  uint32_t radius_buffer = 10;  // Extra radius for RadiusGdtwAligner
};

/// ChainAligner: orchestrates segment alignment for an entire chain.
///
/// Given a chain of anchors and a query signal, aligns each segment
/// between consecutive anchors and aggregates the results.
class ChainAligner {
 public:
  explicit ChainAligner(ChainAlignerConfig config = {});

  /// Align a chain of anchors against the graph.
  ///
  /// @param graph      Graph topology
  /// @param signals    Reference signals per node
  /// @param query      Full query signal (will be sliced per segment)
  /// @param anchors    Ordered anchor points (at least 2)
  /// @return           Aggregated alignment result
  ChainAlignmentResult align(const index::GraphStore& graph,
                             const index::SignalStore& signals,
                             const signal::AlignmentQuantizedSignal& query,
                             const std::vector<Anchor>& anchors);

 private:
  ChainAlignerConfig config_;
  std::unique_ptr<SegmentAligner> path_guided_aligner_;
  std::unique_ptr<SegmentAligner> radius_aligner_;

  /// Select appropriate backend for a segment.
  SegmentAligner* selectBackend(const Anchor& start, const Anchor& target);
};

}  // namespace piru::alignment
