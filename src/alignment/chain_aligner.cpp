// SPDX-License-Identifier: MIT
// ChainAligner implementation.

#include "alignment/chain_aligner.hpp"

#include <limits>

namespace piru::alignment {

ChainAligner::ChainAligner(ChainAlignerConfig config) : config_(config) {
  // Create aligners based on backend configuration
  if (config_.backend == AlignerBackend::kPathGuided ||
      config_.backend == AlignerBackend::kAuto) {
    path_guided_aligner_ = makePathGuidedDtwAligner();
  }

  if (config_.backend == AlignerBackend::kRadius ||
      config_.backend == AlignerBackend::kAuto) {
    radius_aligner_ = makeRadiusGdtwAligner(config_.radius_buffer);
  }
}

SegmentAligner* ChainAligner::selectBackend(const Anchor& start,
                                             const Anchor& target) {
  switch (config_.backend) {
    case AlignerBackend::kPathGuided:
      return path_guided_aligner_.get();

    case AlignerBackend::kRadius:
      return radius_aligner_.get();

    case AlignerBackend::kAuto:
      // Use PathGuided for same-node segments, Radius for cross-node
      if (start.graph_pos.node_id == target.graph_pos.node_id) {
        return path_guided_aligner_.get();
      } else {
        return radius_aligner_.get();
      }
  }
  return path_guided_aligner_.get();  // Fallback
}

ChainAlignmentResult ChainAligner::align(
    const index::GraphStore& graph,
    const index::SignalStore& signals,
    const signal::AlignmentQuantizedSignal& query,
    const std::vector<Anchor>& anchors) {

  ChainAlignmentResult result;
  result.total_cost = 0.0f;
  result.segments_aligned = 0;
  result.segments_reached_target = 0;

  // Need at least 2 anchors to form a segment
  if (anchors.size() < 2) {
    result.total_cost = std::numeric_limits<float>::infinity();
    return result;
  }

  size_t query_len = signalLength(query);

  // Iterate through anchor pairs
  for (size_t i = 0; i + 1 < anchors.size(); ++i) {
    const Anchor& start = anchors[i];
    const Anchor& target = anchors[i + 1];

    // Validate anchor positions
    if (start.query_pos >= query_len || target.query_pos > query_len) {
      continue;  // Skip invalid anchors
    }
    if (start.query_pos >= target.query_pos) {
      continue;  // Skip non-forward anchors
    }

    // Slice query signal for this segment
    // For intermediate segments (not the last), exclude last element
    // to avoid double-counting at anchor boundaries.
    // The anchor position is included in the next segment.
    size_t slice_end = target.query_pos;
    bool is_last_segment = (i + 2 >= anchors.size());
    if (!is_last_segment) {
      // Exclude the anchor position itself for intermediate segments
      // (it will be the start of the next segment)
    }
    // Note: We keep slice_end as target.query_pos for now.
    // The boundary handling is done in segment aligner (target offset is exclusive).

    auto query_segment = sliceSignal(query, start.query_pos, slice_end);

    if (signalLength(query_segment) == 0) {
      continue;  // Skip empty segments
    }

    // Select backend and align
    SegmentAligner* aligner = selectBackend(start, target);
    if (!aligner) {
      continue;  // No aligner available
    }

    DtwResult segment_result = aligner->align(graph, signals, query_segment, start, target);

    if (segment_result.valid()) {
      result.total_cost += segment_result.cost;
      result.segments_aligned++;

      if (segment_result.reached_target) {
        result.segments_reached_target++;
      }

      // Append path, avoiding duplicates at boundaries
      if (!segment_result.path.empty()) {
        if (!result.path.empty() && !is_last_segment) {
          // Skip first position if it duplicates last position from previous segment
          size_t start_idx = 0;
          if (result.path.back() == segment_result.path.front()) {
            start_idx = 1;
          }
          for (size_t j = start_idx; j < segment_result.path.size(); ++j) {
            result.path.push_back(segment_result.path[j]);
          }
        } else {
          // First segment or last segment: append all
          result.path.insert(result.path.end(),
                             segment_result.path.begin(),
                             segment_result.path.end());
        }
      }
    } else {
      // Segment alignment failed - mark chain as partially failed
      // Continue with remaining segments
    }
  }

  // If no segments aligned, mark as invalid
  if (result.segments_aligned == 0) {
    result.total_cost = std::numeric_limits<float>::infinity();
  }

  return result;
}

}  // namespace piru::alignment
