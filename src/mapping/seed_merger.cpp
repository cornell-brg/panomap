/**
 * seed_merger.cpp
 *
 * Merges adjacent seed hits.
 *
 * Related:
 *  - seed_merger.hpp
 *
 * SPDX-License-Identifier: MIT
 */

#include "mapping/seed_merger.hpp"

#include <algorithm>
#include <vector>

namespace piru::mapping {

namespace {

// Comparator for sorting seed hits by (node_id, node_offset, query_pos).
struct SeedHitComparator {
  bool operator()(const NodeAnchor& a, const NodeAnchor& b) const {
    if (a.node_id != b.node_id) return a.node_id < b.node_id;
    if (a.offset != b.offset) return a.offset < b.offset;
    return a.read_pos < b.read_pos;
  }
};

// Check if two seed hits can be merged based on tolerance.
bool can_merge(const NodeAnchor& a, const NodeAnchor& b, std::size_t tolerance) {
  if (a.node_id != b.node_id) return false;

  // Gap between END of a and START of b (sorted order, so b comes after a)
  auto a_query_end = static_cast<std::int64_t>(a.read_pos) + a.span;
  auto a_ref_end = static_cast<std::int64_t>(a.offset) + a.length;

  auto query_gap = static_cast<std::int64_t>(b.read_pos) - a_query_end;
  auto ref_gap = static_cast<std::int64_t>(b.offset) - a_ref_end;

  // Overlap counts as zero gap
  auto qg = (query_gap > 0) ? static_cast<std::size_t>(query_gap) : 0;
  auto rg = (ref_gap > 0) ? static_cast<std::size_t>(ref_gap) : 0;

  return qg <= tolerance && rg <= tolerance;
}

// Merge hit b into hit a (updates a's span to cover both).
void merge_into(NodeAnchor& a, const NodeAnchor& b) {
  auto a_query_end = static_cast<std::uint32_t>(a.read_pos + a.span);
  auto b_query_end = static_cast<std::uint32_t>(b.read_pos + b.span);
  auto a_ref_end = static_cast<std::uint32_t>(a.offset + a.length);
  auto b_ref_end = static_cast<std::uint32_t>(b.offset + b.length);

  std::uint32_t merged_start = std::min(a.read_pos, b.read_pos);
  std::uint32_t merged_end = std::max(a_query_end, b_query_end);
  a.read_pos = merged_start;
  a.span = static_cast<std::uint16_t>(std::min(merged_end - merged_start, 0xFFFFu));

  // Extend ref length to cover both
  std::uint32_t merged_ref_end = std::max(a_ref_end, b_ref_end);
  a.length = static_cast<std::uint16_t>(std::min(merged_ref_end - a.offset, 0xFFFFu));
}

}  // namespace

std::vector<NodeAnchor> SeedMerger::merge(const std::vector<NodeAnchor>& hits,
                                          const SeedMergerConfig& config) {
  if (hits.empty()) {
    return {};
  }

  /* Sort hits by (node_id, node_offset, query_pos) */
  std::vector<NodeAnchor> sorted_hits = hits;
  std::sort(sorted_hits.begin(), sorted_hits.end(), SeedHitComparator{});

  std::vector<NodeAnchor> merged;
  merged.reserve(sorted_hits.size());

  /* Start with the first hit as the current accumulator */
  NodeAnchor current = sorted_hits[0];

  for (std::size_t i = 1; i < sorted_hits.size(); ++i) {
    const auto& next = sorted_hits[i];

    if (can_merge(current, next, config.merge_tolerance)) {
      merge_into(current, next);
    } else {
      merged.push_back(current);
      current = next;
    }
  }

  // last accumulated hit not yet pushed
  merged.push_back(current);

  return merged;
}

}  // namespace piru::mapping
