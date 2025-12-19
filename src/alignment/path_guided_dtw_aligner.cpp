// SPDX-License-Identifier: MIT
// PathGuidedDtwAligner: linear extraction + standard DTW.

#include "alignment/segment_aligner.hpp"
#include "alignment/signal_utils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <variant>

namespace piru::alignment {
namespace {

// =============================================================================
// Linear Extraction
// =============================================================================

/// Find path of nodes from start_node to end_node using BFS.
/// Returns empty vector if no path found.
std::vector<size_t> findNodePath(const index::GraphStore& graph,
                                  size_t start_node, size_t end_node) {
  if (start_node == end_node) {
    return {start_node};
  }

  // BFS to find path
  std::queue<size_t> queue;
  std::unordered_map<size_t, size_t> parent;  // child -> parent

  queue.push(start_node);
  parent[start_node] = start_node;  // sentinel

  while (!queue.empty()) {
    size_t curr = queue.front();
    queue.pop();

    if (curr == end_node) {
      // Reconstruct path
      std::vector<size_t> path;
      size_t node = end_node;
      while (node != start_node) {
        path.push_back(node);
        node = parent[node];
      }
      path.push_back(start_node);
      std::reverse(path.begin(), path.end());
      return path;
    }

    for (size_t next : graph.outgoing(curr)) {
      if (parent.find(next) == parent.end()) {
        parent[next] = curr;
        queue.push(next);
      }
    }
  }

  return {};  // No path found
}

/// Extract linear signal along path between anchors.
/// Returns concatenated signal from start offset to target offset.
signal::AlignmentQuantizedSignal extractLinearSignal(
    const index::GraphStore& graph,
    const index::SignalStore& signals,
    const Anchor& start,
    const Anchor& target) {

  std::vector<size_t> node_path = findNodePath(graph, start.graph_pos.node_id,
                                                target.graph_pos.node_id);

  if (node_path.empty()) {
    // No path found - return empty signal
    signal::AlignmentQuantizedSignal empty;
    empty.kind = signal::AlignmentQuantizationKind::kFloat32;
    empty.data = std::vector<float>{};
    return empty;
  }

  // Get first node's signal to determine quantization kind
  const auto* first_sig = signals.get(node_path[0]);
  if (!first_sig) {
    signal::AlignmentQuantizedSignal empty;
    empty.kind = signal::AlignmentQuantizationKind::kFloat32;
    empty.data = std::vector<float>{};
    return empty;
  }

  signal::AlignmentQuantizedSignal result;
  result.kind = first_sig->kind;

  // Initialize with empty vector of correct type
  std::visit([&result](const auto& vec) {
    using T = typename std::decay_t<decltype(vec)>::value_type;
    result.data = std::vector<T>{};
  }, first_sig->data);

  for (size_t i = 0; i < node_path.size(); ++i) {
    size_t node_id = node_path[i];
    const auto* node_sig = signals.get(node_id);
    if (!node_sig) continue;

    size_t sig_len = signalLength(*node_sig);
    size_t slice_start = 0;
    size_t slice_end = sig_len;

    if (i == 0) {
      // First node: from start offset
      slice_start = start.graph_pos.offset;
    }
    if (i == node_path.size() - 1) {
      // Last node: to target offset (exclusive for boundary handling)
      slice_end = target.graph_pos.offset;
    }

    if (slice_start < slice_end) {
      auto sliced = sliceSignal(*node_sig, slice_start, slice_end);
      result = concatSignals(result, sliced);
    }
  }

  return result;
}

// =============================================================================
// Linear DTW
// =============================================================================

/// Compute absolute distance between two signal values.
template <typename T>
inline float distance(T a, T b) {
  return std::abs(static_cast<float>(a) - static_cast<float>(b));
}

/// Standard linear DTW with traceback.
/// Returns cost and path (indices into query and ref).
template <typename T>
std::pair<float, std::vector<std::pair<size_t, size_t>>> linearDtw(
    const std::vector<T>& query,
    const std::vector<T>& ref) {

  if (query.empty() || ref.empty()) {
    return {std::numeric_limits<float>::infinity(), {}};
  }

  size_t m = query.size();
  size_t n = ref.size();

  // DP matrix
  std::vector<std::vector<float>> dp(m, std::vector<float>(n, std::numeric_limits<float>::infinity()));

  // Parent pointers: 0 = diagonal, 1 = vertical (from i-1), 2 = horizontal (from j-1)
  std::vector<std::vector<uint8_t>> parent(m, std::vector<uint8_t>(n, 0));

  // Initialize
  dp[0][0] = distance(query[0], ref[0]);

  // First column (vertical moves only)
  for (size_t i = 1; i < m; ++i) {
    dp[i][0] = dp[i-1][0] + distance(query[i], ref[0]);
    parent[i][0] = 1;  // came from above
  }

  // First row (horizontal moves only)
  for (size_t j = 1; j < n; ++j) {
    dp[0][j] = dp[0][j-1] + distance(query[0], ref[j]);
    parent[0][j] = 2;  // came from left
  }

  // Fill DP matrix
  for (size_t i = 1; i < m; ++i) {
    for (size_t j = 1; j < n; ++j) {
      float cost = distance(query[i], ref[j]);
      float diag = dp[i-1][j-1];
      float vert = dp[i-1][j];
      float horz = dp[i][j-1];

      if (diag <= vert && diag <= horz) {
        dp[i][j] = diag + cost;
        parent[i][j] = 0;
      } else if (vert <= horz) {
        dp[i][j] = vert + cost;
        parent[i][j] = 1;
      } else {
        dp[i][j] = horz + cost;
        parent[i][j] = 2;
      }
    }
  }

  // Traceback
  std::vector<std::pair<size_t, size_t>> path;
  size_t i = m - 1;
  size_t j = n - 1;

  while (i > 0 || j > 0) {
    path.push_back({i, j});
    uint8_t p = parent[i][j];
    if (p == 0) {
      --i; --j;
    } else if (p == 1) {
      --i;
    } else {
      --j;
    }
  }
  path.push_back({0, 0});
  std::reverse(path.begin(), path.end());

  return {dp[m-1][n-1], path};
}

/// Dispatch DTW based on signal type.
std::pair<float, std::vector<std::pair<size_t, size_t>>> runDtw(
    const signal::AlignmentQuantizedSignal& query,
    const signal::AlignmentQuantizedSignal& ref) {

  if (query.kind != ref.kind) {
    return {std::numeric_limits<float>::infinity(), {}};
  }

  return std::visit([&ref](const auto& query_vec) {
    using T = typename std::decay_t<decltype(query_vec)>::value_type;
    const auto& ref_vec = std::get<std::vector<T>>(ref.data);
    return linearDtw(query_vec, ref_vec);
  }, query.data);
}

// =============================================================================
// PathGuidedDtwAligner Implementation
// =============================================================================

class PathGuidedDtwAligner : public SegmentAligner {
 public:
  DtwResult align(const index::GraphStore& graph,
                  const index::SignalStore& signals,
                  const signal::AlignmentQuantizedSignal& query,
                  Anchor start, Anchor target) override {

    // Extract reference signal along path
    auto ref_signal = extractLinearSignal(graph, signals, start, target);

    if (signalLength(ref_signal) == 0 || signalLength(query) == 0) {
      return DtwResult{
          .cost = std::numeric_limits<float>::infinity(),
          .end_pos = start.graph_pos,
          .reached_target = false,
          .path = {}
      };
    }

    // Run DTW
    auto [cost, dtw_path] = runDtw(query, ref_signal);

    // Convert DTW path (ref indices) to graph positions
    // For linear extraction, we need to map ref index back to graph position
    std::vector<GraphPosition> graph_path;
    graph_path.reserve(dtw_path.size());

    // Build mapping from linear ref index to graph position
    std::vector<GraphPosition> ref_index_to_pos;
    {
      auto node_path = findNodePath(graph, start.graph_pos.node_id,
                                    target.graph_pos.node_id);

      for (size_t i = 0; i < node_path.size(); ++i) {
        size_t node_id = node_path[i];
        const auto* node_sig = signals.get(node_id);
        if (!node_sig) continue;

        size_t sig_len = signalLength(*node_sig);
        size_t slice_start = 0;
        size_t slice_end = sig_len;

        if (i == 0) slice_start = start.graph_pos.offset;
        if (i == node_path.size() - 1) slice_end = target.graph_pos.offset;

        for (size_t off = slice_start; off < slice_end; ++off) {
          ref_index_to_pos.push_back(GraphPosition{
              static_cast<uint32_t>(node_id),
              static_cast<uint32_t>(off)
          });
        }
      }
    }

    // Map DTW path to graph positions
    for (const auto& [q_idx, r_idx] : dtw_path) {
      if (r_idx < ref_index_to_pos.size()) {
        graph_path.push_back(ref_index_to_pos[r_idx]);
      }
    }

    GraphPosition end_pos = graph_path.empty() ? start.graph_pos : graph_path.back();
    bool reached = (end_pos.node_id == target.graph_pos.node_id &&
                    end_pos.offset + 1 >= target.graph_pos.offset);  // +1 for boundary

    return DtwResult{
        .cost = cost,
        .end_pos = end_pos,
        .reached_target = reached,
        .path = std::move(graph_path)
    };
  }

  std::string name() const override { return "PathGuidedDtwAligner"; }
};

}  // namespace

// Factory function
std::unique_ptr<SegmentAligner> makePathGuidedDtwAligner() {
  return std::make_unique<PathGuidedDtwAligner>();
}

}  // namespace piru::alignment
