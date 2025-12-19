// SPDX-License-Identifier: MIT
// RadiusGdtwAligner: BFS extraction + Navarro-style graph DTW.

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
// Window: Extracted subgraph for alignment
// =============================================================================

struct Window {
  // Positions in the window (local ID -> graph position)
  std::vector<GraphPosition> positions;

  // Signals at each position (local ID -> signal value as float)
  std::vector<float> signals;

  // Edges: successors[k] contains local IDs of positions reachable from k
  std::vector<std::vector<uint32_t>> successors;

  // Special indices
  uint32_t start_idx = 0;
  int32_t target_idx = -1;  // -1 if target not in window

  size_t size() const { return positions.size(); }
};

// =============================================================================
// BFS Radius Extraction
// =============================================================================

/// Extract window via BFS from start position.
/// Includes all positions within radius steps.
Window extractWindow(const index::GraphStore& graph,
                     const index::SignalStore& signals,
                     const Anchor& start,
                     const Anchor& target,
                     uint32_t radius) {
  Window window;

  // Map from (node_id, offset) to local ID
  std::unordered_map<uint64_t, uint32_t> visited;
  auto posKey = [](uint32_t node, uint32_t off) -> uint64_t {
    return (static_cast<uint64_t>(node) << 32) | off;
  };

  // BFS queue: (node_id, offset, distance)
  std::queue<std::tuple<uint32_t, uint32_t, uint32_t>> queue;
  queue.push({start.graph_pos.node_id, start.graph_pos.offset, 0});

  while (!queue.empty()) {
    auto [node_id, offset, dist] = queue.front();
    queue.pop();

    if (dist > radius) continue;

    uint64_t key = posKey(node_id, offset);
    if (visited.count(key)) continue;

    // Get signal for this node
    const auto* node_sig = signals.get(node_id);
    if (!node_sig) continue;

    size_t sig_len = signalLength(*node_sig);
    if (offset >= sig_len) continue;

    // Add to window
    uint32_t local_id = static_cast<uint32_t>(window.positions.size());
    visited[key] = local_id;
    window.positions.push_back(GraphPosition{node_id, offset});
    window.signals.push_back(signalValueAt(*node_sig, offset));
    window.successors.push_back({});

    // Check if this is start or target
    if (node_id == start.graph_pos.node_id && offset == start.graph_pos.offset) {
      window.start_idx = local_id;
    }
    if (node_id == target.graph_pos.node_id && offset == target.graph_pos.offset) {
      window.target_idx = static_cast<int32_t>(local_id);
    }

    // Enqueue successors
    if (offset + 1 < sig_len) {
      // Next position within same node
      queue.push({node_id, offset + 1, dist + 1});
    } else {
      // Transition to successor nodes (offset 0)
      for (size_t succ_id : graph.outgoing(node_id)) {
        queue.push({static_cast<uint32_t>(succ_id), 0, dist + 1});
      }
    }
  }

  // Build edge lists (successors for each position)
  for (uint32_t k = 0; k < window.positions.size(); ++k) {
    const auto& pos = window.positions[k];
    const auto* node_sig = signals.get(pos.node_id);
    if (!node_sig) continue;

    size_t sig_len = signalLength(*node_sig);

    if (pos.offset + 1 < sig_len) {
      // Successor within same node
      uint64_t succ_key = posKey(pos.node_id, pos.offset + 1);
      auto it = visited.find(succ_key);
      if (it != visited.end()) {
        window.successors[k].push_back(it->second);
      }
    } else {
      // Successors in next nodes
      for (size_t succ_node : graph.outgoing(pos.node_id)) {
        uint64_t succ_key = posKey(static_cast<uint32_t>(succ_node), 0);
        auto it = visited.find(succ_key);
        if (it != visited.end()) {
          window.successors[k].push_back(it->second);
        }
      }
    }
  }

  return window;
}

// =============================================================================
// Navarro-style Graph DTW
// =============================================================================

/// Compute distance between query value and reference value.
inline float dist(float q, float r) {
  return std::abs(q - r);
}

/// Two-phase graph DTW: vertical+diagonal, then horizontal propagation.
/// Returns (cost, path as local IDs).
std::pair<float, std::vector<uint32_t>> graphDtw(
    const Window& window,
    const std::vector<float>& query) {

  if (query.empty() || window.size() == 0) {
    return {std::numeric_limits<float>::infinity(), {}};
  }

  size_t m = query.size();
  size_t K = window.size();
  constexpr float INF = std::numeric_limits<float>::infinity();

  // DP arrays: D_prev[k] = cost to reach position k after processing query[0..i-1]
  std::vector<float> D_prev(K, INF);
  std::vector<float> D_cur(K, INF);

  // Parent tracking for traceback: parent[i][k] = (prev_query_idx, prev_local_id)
  // We store parents per row for memory efficiency
  std::vector<std::vector<std::pair<int32_t, int32_t>>> parent(m, std::vector<std::pair<int32_t, int32_t>>(K, {-1, -1}));

  // Initialize: only start position is reachable at row 0
  D_prev[window.start_idx] = dist(query[0], window.signals[window.start_idx]);

  // Horizontal propagation for row 0
  {
    std::queue<uint32_t> prop_queue;
    prop_queue.push(window.start_idx);
    std::vector<bool> in_queue(K, false);
    in_queue[window.start_idx] = true;

    while (!prop_queue.empty()) {
      uint32_t x = prop_queue.front();
      prop_queue.pop();
      in_queue[x] = false;

      for (uint32_t y : window.successors[x]) {
        float emit = dist(query[0], window.signals[y]);
        float candidate = D_prev[x] + emit;
        if (candidate < D_prev[y]) {
          D_prev[y] = candidate;
          parent[0][y] = {0, static_cast<int32_t>(x)};
          if (!in_queue[y]) {
            prop_queue.push(y);
            in_queue[y] = true;
          }
        }
      }
    }
  }

  // Main DP loop
  for (size_t i = 1; i < m; ++i) {
    std::fill(D_cur.begin(), D_cur.end(), INF);

    // Phase 1: Vertical + Diagonal (no same-row dependencies)
    for (uint32_t k = 0; k < K; ++k) {
      if (D_prev[k] >= INF) continue;

      float emit_k = dist(query[i], window.signals[k]);

      // Vertical: stay at position k, advance query
      if (D_prev[k] + emit_k < D_cur[k]) {
        D_cur[k] = D_prev[k] + emit_k;
        parent[i][k] = {static_cast<int32_t>(i - 1), static_cast<int32_t>(k)};
      }

      // Diagonal: advance to successor k', advance query
      for (uint32_t kp : window.successors[k]) {
        float emit_kp = dist(query[i], window.signals[kp]);
        if (D_prev[k] + emit_kp < D_cur[kp]) {
          D_cur[kp] = D_prev[k] + emit_kp;
          parent[i][kp] = {static_cast<int32_t>(i - 1), static_cast<int32_t>(k)};
        }
      }
    }

    // Phase 2: Horizontal propagation (handles cycles)
    {
      std::queue<uint32_t> prop_queue;
      std::vector<bool> in_queue(K, false);

      for (uint32_t k = 0; k < K; ++k) {
        if (D_cur[k] < INF) {
          prop_queue.push(k);
          in_queue[k] = true;
        }
      }

      while (!prop_queue.empty()) {
        uint32_t x = prop_queue.front();
        prop_queue.pop();
        in_queue[x] = false;

        for (uint32_t y : window.successors[x]) {
          float emit_y = dist(query[i], window.signals[y]);
          float candidate = D_cur[x] + emit_y;
          if (candidate < D_cur[y]) {
            D_cur[y] = candidate;
            parent[i][y] = {static_cast<int32_t>(i), static_cast<int32_t>(x)};
            if (!in_queue[y]) {
              prop_queue.push(y);
              in_queue[y] = true;
            }
          }
        }
      }
    }

    std::swap(D_prev, D_cur);
  }

  // Find best endpoint
  uint32_t best_k = 0;
  float best_cost = INF;
  for (uint32_t k = 0; k < K; ++k) {
    if (D_prev[k] < best_cost) {
      best_cost = D_prev[k];
      best_k = k;
    }
  }

  if (best_cost >= INF) {
    return {INF, {}};
  }

  // Traceback
  std::vector<uint32_t> path;
  int32_t qi = static_cast<int32_t>(m - 1);
  int32_t ki = static_cast<int32_t>(best_k);

  while (qi >= 0 && ki >= 0) {
    path.push_back(static_cast<uint32_t>(ki));
    auto [prev_q, prev_k] = parent[qi][ki];
    qi = prev_q;
    ki = prev_k;
  }

  std::reverse(path.begin(), path.end());

  return {best_cost, path};
}

// =============================================================================
// RadiusGdtwAligner Implementation
// =============================================================================

class RadiusGdtwAligner : public SegmentAligner {
 public:
  explicit RadiusGdtwAligner(uint32_t radius_buffer) : radius_buffer_(radius_buffer) {}

  DtwResult align(const index::GraphStore& graph,
                  const index::SignalStore& signals,
                  const signal::AlignmentQuantizedSignal& query,
                  Anchor start, Anchor target) override {

    // Convert query to float vector
    std::vector<float> query_floats;
    std::visit([&query_floats](const auto& vec) {
      query_floats.reserve(vec.size());
      for (const auto& v : vec) {
        query_floats.push_back(static_cast<float>(v));
      }
    }, query.data);

    if (query_floats.empty()) {
      return DtwResult{
          .cost = std::numeric_limits<float>::infinity(),
          .end_pos = start.graph_pos,
          .reached_target = false,
          .path = {}
      };
    }

    // Compute radius: expected distance + buffer
    // Expected distance approximated by query length (assuming ~1:1 rate)
    uint32_t radius = static_cast<uint32_t>(query_floats.size()) + radius_buffer_;

    // Extract window via BFS
    Window window = extractWindow(graph, signals, start, target, radius);

    if (window.size() == 0) {
      return DtwResult{
          .cost = std::numeric_limits<float>::infinity(),
          .end_pos = start.graph_pos,
          .reached_target = false,
          .path = {}
      };
    }

    // Run graph DTW
    auto [cost, local_path] = graphDtw(window, query_floats);

    // Convert local path to graph positions
    std::vector<GraphPosition> graph_path;
    graph_path.reserve(local_path.size());
    for (uint32_t local_id : local_path) {
      if (local_id < window.positions.size()) {
        graph_path.push_back(window.positions[local_id]);
      }
    }

    GraphPosition end_pos = graph_path.empty() ? start.graph_pos : graph_path.back();

    // Check if we reached target
    bool reached = (window.target_idx >= 0) &&
                   !local_path.empty() &&
                   (local_path.back() == static_cast<uint32_t>(window.target_idx));

    return DtwResult{
        .cost = cost,
        .end_pos = end_pos,
        .reached_target = reached,
        .path = std::move(graph_path)
    };
  }

  std::string name() const override { return "RadiusGdtwAligner"; }

 private:
  uint32_t radius_buffer_;
};

}  // namespace

// Factory function
std::unique_ptr<SegmentAligner> makeRadiusGdtwAligner(uint32_t radius_buffer) {
  return std::make_unique<RadiusGdtwAligner>(radius_buffer);
}

}  // namespace piru::alignment
