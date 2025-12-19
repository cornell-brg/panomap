// SPDX-License-Identifier: MIT
// Segment alignment interface for chain evaluation.

#pragma once

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "index/graph_store.hpp"
#include "index/signal_store.hpp"
#include "signal/signal_types.hpp"

namespace piru::alignment {

/// Position in the graph (node + offset within node's signal).
struct GraphPosition {
  uint32_t node_id;
  uint32_t offset;  // Position within node's signal array

  bool operator==(const GraphPosition& other) const {
    return node_id == other.node_id && offset == other.offset;
  }
};

/// Anchor point for segment alignment.
struct Anchor {
  GraphPosition graph_pos;
  uint32_t query_pos;  // Position in query signal
};

/// Result of segment alignment.
struct DtwResult {
  float cost;                           // DTW cost (lower = better)
  GraphPosition end_pos;                // Where alignment ended
  bool reached_target;                  // Did we reach target anchor?
  std::vector<GraphPosition> path;      // Alignment path through graph

  bool valid() const { return cost < std::numeric_limits<float>::infinity(); }
};

/// Abstract interface for segment alignment.
///
/// Aligns a query signal segment against reference signals in the graph
/// between two anchor positions.
class SegmentAligner {
 public:
  virtual ~SegmentAligner() = default;

  /// Align query segment between start and target anchors.
  ///
  /// @param graph      Graph topology
  /// @param signals    Reference signals per node
  /// @param query      Query signal segment (already sliced for this segment)
  /// @param start      Start anchor position
  /// @param target     Target anchor position (soft constraint)
  /// @return           Alignment result with cost, end position, and path
  virtual DtwResult align(const index::GraphStore& graph,
                          const index::SignalStore& signals,
                          const signal::AlignmentQuantizedSignal& query,
                          Anchor start, Anchor target) = 0;

  /// Backend name for logging/debugging.
  virtual std::string name() const = 0;
};

// Factory functions (implemented in respective backend source files)

/// Create PathGuidedDtwAligner: linear extraction + standard DTW.
/// For same-path anchors (current chaining constraint).
std::unique_ptr<SegmentAligner> makePathGuidedDtwAligner();

/// Create RadiusGdtwAligner: BFS extraction + Navarro-style graph DTW.
/// For bubble exploration and future cross-path chaining.
/// @param radius_buffer Extra positions beyond expected distance
std::unique_ptr<SegmentAligner> makeRadiusGdtwAligner(uint32_t radius_buffer = 0);

}  // namespace piru::alignment
