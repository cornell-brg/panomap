// SPDX-License-Identifier: MIT

#include "mapping/seed_merger.hpp"

#include <algorithm>
#include <vector>

namespace piru::mapping {

namespace {

// Comparator for sorting seed hits by (node_id, node_offset, query_pos).
struct SeedHitComparator {
    bool operator()(const NodeAnchor& a, const NodeAnchor& b) const {
        if (a.target.node_id != b.target.node_id) {
            return a.target.node_id < b.target.node_id;
        }
        if (a.target.offset != b.target.offset) {
            return a.target.offset < b.target.offset;
        }
        return a.read_pos < b.read_pos;
    }
};

// Check if two seed hits can be merged based on tolerance.
bool can_merge(const NodeAnchor& a, const NodeAnchor& b, std::size_t tolerance) {
    // Must be on same node
    if (a.target.node_id != b.target.node_id) {
        return false;
    }

    // Compute gap between END of a and START of b
    // If b starts before a ends, gap is 0 (overlap)
    // Note: hits are sorted, so b comes after a in sorted order
    const std::size_t a_query_end = a.read_pos + a.span;
    const std::size_t a_ref_end = a.target.offset + a.target.length;

    const std::int64_t query_gap_signed =
        static_cast<std::int64_t>(b.read_pos) - static_cast<std::int64_t>(a_query_end);
    const std::int64_t ref_gap_signed =
        static_cast<std::int64_t>(b.target.offset) - static_cast<std::int64_t>(a_ref_end);

    // Gap is 0 if there's overlap (negative gap)
    const std::size_t query_gap =
        (query_gap_signed > 0) ? static_cast<std::size_t>(query_gap_signed) : 0;
    const std::size_t ref_gap = (ref_gap_signed > 0) ? static_cast<std::size_t>(ref_gap_signed) : 0;

    // Allow merging if gaps are within tolerance
    return query_gap <= tolerance && ref_gap <= tolerance;
}

// Merge hit b into hit a (updates a's span to cover both).
void merge_into(NodeAnchor& a, const NodeAnchor& b) {
    // Calculate the end positions
    const std::size_t a_end = a.read_pos + a.span;
    const std::size_t b_end = b.read_pos + b.span;

    // Extend span to cover both seeds
    const std::size_t merged_start = std::min(a.read_pos, b.read_pos);
    const std::size_t merged_end = std::max(a_end, b_end);

    // Update a to represent the merged hit
    a.read_pos = merged_start;
    a.span = merged_end - merged_start;

    // Use the earlier node offset (from a, which is already sorted first)
    // Keep other fields from a (hash, frequency, score)
}

}  // namespace

std::vector<NodeAnchor> SeedMerger::merge(const std::vector<NodeAnchor>& hits,
                                             const SeedMergerConfig& config) {
    if (hits.empty()) {
        return {};
    }

    // Sort hits by (node_id, node_offset, query_pos)
    std::vector<NodeAnchor> sorted_hits = hits;
    std::sort(sorted_hits.begin(), sorted_hits.end(), SeedHitComparator{});

    std::vector<NodeAnchor> merged;
    merged.reserve(sorted_hits.size());

    // Start with the first hit as the current accumulator
    NodeAnchor current = sorted_hits[0];

    for (std::size_t i = 1; i < sorted_hits.size(); ++i) {
        const auto& next = sorted_hits[i];

        if (can_merge(current, next, config.merge_tolerance)) {
            // Merge next into current
            merge_into(current, next);
        } else {
            // Cannot merge - push current and start new accumulator
            merged.push_back(current);
            current = next;
        }
    }

    // Don't forget the last accumulated hit
    merged.push_back(current);

    return merged;
}

}  // namespace piru::mapping
