// SPDX-License-Identifier: MIT
// Superbubble-based linearization implementation.

#include "index/superbubble_linearizer.hpp"

#include "index/pseudo_linearize.hpp"

namespace piru::index {

std::vector<std::vector<LinearCoordinate>> SuperbubbleLinearizer::linearize(
    const AlnGraph& graph) const {
  const std::size_t n = graph.nodeCount();

  // Run existing pseudo-linearization pipeline.
  SccResult scc = computeScc(graph);
  TipFoldingResult tips = chainTips(graph, scc);
  chainCycles(graph, tips);
  SuperbubbleResult sb = chainSuperbubbles(graph, scc, tips);
  auto chain_ids = assignChainIds(std::move(sb.uf));
  auto positions = assignLinearPositions(graph, chain_ids, scc);

  // Convert to LinearCoordinate format.
  // Each node gets at most one coordinate (its chain_id and linear_position).
  std::vector<std::vector<LinearCoordinate>> result(n);

  for (std::size_t i = 0; i < n; ++i) {
    // Check if node has valid linear position (not -1).
    // Note: assignLinearPositions may return -1 for unassigned nodes.
    if (positions[i] != -1) {
      result[i].emplace_back(chain_ids[i], positions[i]);
    }
    // Otherwise, result[i] remains empty (no coordinates).
  }

  return result;
}

}  // namespace piru::index
