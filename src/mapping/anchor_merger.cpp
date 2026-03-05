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

// Check if two anchors can be merged.
// Rules:
// 1. Must be on same diagonal (exact matches have equal query/ref spans)
// 2. Must overlap (no gap between them)
bool can_merge(const Anchor& a, const Anchor& b) {
    // Must be on same diagonal (diagonal = ref_coord - query_pos)
    const std::int64_t a_diag = a.ref_coord - static_cast<std::int64_t>(a.query_pos);
    const std::int64_t b_diag = b.ref_coord - static_cast<std::int64_t>(b.query_pos);
    if (a_diag != b_diag) {
        return false;  // Different diagonals = different alignments, don't merge
    }

    // On same diagonal, just check if they overlap in query (ref will match)
    const std::size_t a_query_end = a.query_pos + a.length;
    return b.query_pos <= a_query_end;
}

// Merge anchor b into anchor a (updates a's length to cover both).
// Precondition: a and b are on the same diagonal, so query_span == ref_span.
void merge_into(Anchor& a, const Anchor& b) {
    // Since they're on the same diagonal, we only need to track query positions
    // (ref positions will be consistent due to same diagonal)
    const std::size_t a_query_end = a.query_pos + a.length;
    const std::size_t b_query_end = b.query_pos + b.length;

    // Extend to cover both (a starts first since sorted by ref_coord)
    const std::size_t merged_query_end = std::max(a_query_end, b_query_end);
    a.length = merged_query_end - a.query_pos;

    // Keep a's query_pos, ref_coord, path_id, node_id, node_offset (first anchor)
}

}  // namespace

std::vector<Anchor> AnchorMerger::merge(const std::vector<Anchor>& anchors,
                                        const AnchorMergerConfig& /*config*/) {
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

            if (can_merge(current, next)) {
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
