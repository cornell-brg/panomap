// SPDX-License-Identifier: MIT

#include "mapping/anchor_merger.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace piru::mapping {

namespace {

// Comparator for sorting anchors by (ref_coord, query_pos).
struct AnchorComparator {
    bool operator()(const Anchor& a, const Anchor& b) const {
        if (a.ref_coord != b.ref_coord) {
            return a.ref_coord < b.ref_coord;
        }
        return a.query_pos < b.query_pos;
    }
};

// Check if two anchors can be merged based on tolerance.
bool can_merge(const Anchor& a, const Anchor& b, std::size_t tolerance) {
    // Compute gap between END of a and START of b
    const std::int64_t a_query_end = static_cast<std::int64_t>(a.query_pos + a.length);
    const std::int64_t a_ref_end = a.ref_coord + static_cast<std::int64_t>(a.length);

    const std::int64_t query_gap_signed = static_cast<std::int64_t>(b.query_pos) - a_query_end;
    const std::int64_t ref_gap_signed = b.ref_coord - a_ref_end;

    // Reject merging if query position goes backwards (large negative gap)
    // Allow small overlaps (within tolerance) but not large backward jumps
    if (query_gap_signed < -static_cast<std::int64_t>(tolerance)) {
        return false;
    }
    if (ref_gap_signed < -static_cast<std::int64_t>(tolerance)) {
        return false;
    }

    // Compute absolute gap (0 for overlaps within tolerance)
    const std::size_t query_gap = (query_gap_signed > 0) ?
                                  static_cast<std::size_t>(query_gap_signed) : 0;
    const std::size_t ref_gap = (ref_gap_signed > 0) ?
                                static_cast<std::size_t>(ref_gap_signed) : 0;

    // Allow merging if forward gaps are within tolerance
    return query_gap <= tolerance && ref_gap <= tolerance;
}

// Merge anchor b into anchor a (updates a's length to cover both).
void merge_into(Anchor& a, const Anchor& b) {
    // Calculate the end positions
    const std::size_t a_query_end = a.query_pos + a.length;
    const std::size_t b_query_end = b.query_pos + b.length;
    const std::int64_t a_ref_end = a.ref_coord + static_cast<std::int64_t>(a.length);
    const std::int64_t b_ref_end = b.ref_coord + static_cast<std::int64_t>(b.length);

    // Extend span to cover both anchors
    const std::size_t merged_query_start = std::min(a.query_pos, b.query_pos);
    const std::size_t merged_query_end = std::max(a_query_end, b_query_end);
    const std::int64_t merged_ref_start = std::min(a.ref_coord, b.ref_coord);
    const std::int64_t merged_ref_end = std::max(a_ref_end, b_ref_end);

    // Update anchor to represent merged span
    // Use the maximum span needed on either query or reference dimension
    a.query_pos = merged_query_start;
    a.ref_coord = merged_ref_start;
    a.length = std::max(merged_query_end - merged_query_start,
                        static_cast<std::size_t>(merged_ref_end - merged_ref_start));

    // Keep other fields from a (path_id, node_id, node_offset from first anchor)
}

}  // namespace

std::vector<Anchor> AnchorMerger::merge(
    const std::vector<Anchor>& anchors,
    const AnchorMergerConfig& config) {

    if (anchors.empty()) {
        return {};
    }

    // Group anchors by path_id
    std::map<std::size_t, std::vector<Anchor>> anchors_by_path;
    for (const auto& anchor : anchors) {
        anchors_by_path[anchor.path_id].push_back(anchor);
    }

    std::vector<Anchor> merged;
    merged.reserve(anchors.size());  // Upper bound on output size

    // Process each path independently
    for (auto& [path_id, path_anchors] : anchors_by_path) {
        // Sort by (ref_coord, query_pos) within this path
        std::sort(path_anchors.begin(), path_anchors.end(), AnchorComparator{});

        // Merge consecutive anchors
        Anchor current = path_anchors[0];

        for (std::size_t i = 1; i < path_anchors.size(); ++i) {
            const auto& next = path_anchors[i];

            if (can_merge(current, next, config.merge_tolerance)) {
                // Merge next into current
                merge_into(current, next);
            } else {
                // Cannot merge - push current and start new accumulator
                merged.push_back(current);
                current = next;
            }
        }

        // Don't forget the last accumulated anchor
        merged.push_back(current);
    }

    return merged;
}

}  // namespace piru::mapping
