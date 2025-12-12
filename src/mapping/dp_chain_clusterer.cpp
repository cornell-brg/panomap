// SPDX-License-Identifier: MIT

#include "mapping/dp_chain_clusterer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace piru::mapping {

namespace {

// Comparator for sorting anchors by (path_id, ref_coord, query_pos).
struct AnchorComparator {
    bool operator()(const Anchor& a, const Anchor& b) const {
        if (a.path_id != b.path_id) {
            return a.path_id < b.path_id;
        }
        if (a.ref_coord != b.ref_coord) {
            return a.ref_coord < b.ref_coord;
        }
        return a.query_pos < b.query_pos;
    }
};

}  // namespace

DPChainClusterer::DPChainClusterer(
    const std::vector<std::vector<index::LinearCoordinate>>& coords,
    DPChainClustererConfig config)
    : expander_(coords), config_(std::move(config)) {}

ClusterSummary DPChainClusterer::cluster(const std::vector<SeedHitRecord>& hits) const {
    ClusterSummary summary;

    if (hits.empty()) {
        return summary;
    }

    // Step 1: Expand seed hits to anchors
    auto anchors = expander_.expand(hits);

    if (anchors.empty()) {
        return summary;
    }

    // Step 2: Sort anchors by (path_id, ref_coord, query_pos)
    std::sort(anchors.begin(), anchors.end(), AnchorComparator{});

    const std::size_t n = anchors.size();

    // Step 3: DP initialization
    std::vector<double> dp(n, 0.0);
    std::vector<int> pred(n, -1);  // -1 means no predecessor (chain starts here)

    // Step 4: DP loop - compute best score for each anchor
    for (std::size_t i = 0; i < n; ++i) {
        const auto& anchor_i = anchors[i];
        double best_score = anchor_score(anchor_i);  // Start new chain at i
        int best_pred = -1;

        // Check all potential predecessors j < i
        for (std::size_t j = 0; j < i; ++j) {
            const auto& anchor_j = anchors[j];

            // Check if j can chain to i
            if (!can_chain(anchor_j, anchor_i)) {
                continue;
            }

            // Compute score if we extend chain ending at j with anchor i
            const double cost = gap_cost(anchor_j, anchor_i);
            const double score = dp[j] + anchor_score(anchor_i) - cost;

            if (score > best_score) {
                best_score = score;
                best_pred = static_cast<int>(j);
            }
        }

        dp[i] = best_score;
        pred[i] = best_pred;
    }

    // Step 5: Find anchor with best DP score
    std::size_t best_idx = 0;
    double best_dp_score = dp[0];
    for (std::size_t i = 1; i < n; ++i) {
        if (dp[i] > best_dp_score) {
            best_dp_score = dp[i];
            best_idx = i;
        }
    }

    // Check if best chain meets minimum score threshold
    if (best_dp_score < static_cast<double>(config_.min_chain_score)) {
        return summary;  // No chain meets threshold
    }

    // Step 6: Backtrack to extract chain
    auto chain_indices = backtrack_chain(pred, best_idx);

    // Step 7: Convert chain anchors to SeedAnchor format
    summary.score = best_dp_score;
    summary.anchors.reserve(chain_indices.size());

    for (std::size_t idx : chain_indices) {
        const auto& anchor = anchors[idx];
        SeedAnchor seed_anchor;
        seed_anchor.target.node_id = anchor.node_id;
        seed_anchor.target.offset = anchor.node_offset;
        seed_anchor.target.length = anchor.length;
        seed_anchor.read_pos = anchor.query_pos;
        seed_anchor.score = dp[idx];
        seed_anchor.cluster_id = 0;  // Single chain for now

        summary.anchors.push_back(seed_anchor);
    }

    // TODO: Multi-chain extraction (Phase 5 optional subtask)
    // For now, just return the best chain

    return summary;
}

bool DPChainClusterer::can_chain(const Anchor& j, const Anchor& i) const {
    // Check order: j must come before i in sorted order (already guaranteed by DP loop)
    // But also check that query positions are in order
    if (i.query_pos < j.query_pos) {
        return false;
    }

    // Check path constraint
    if (!config_.allow_cross_haplotypes && i.path_id != j.path_id) {
        return false;
    }

    // Compute deltas
    const std::int64_t delta_ref = i.ref_coord - j.ref_coord;
    const std::int64_t delta_query = static_cast<std::int64_t>(i.query_pos) -
                                     static_cast<std::int64_t>(j.query_pos);

    // Check distance constraints
    if (delta_ref < 0 || delta_query < 0) {
        return false;  // Must be forward in both dimensions
    }

    if (static_cast<std::size_t>(delta_ref) > config_.max_dist ||
        static_cast<std::size_t>(delta_query) > config_.max_dist) {
        return false;
    }

    // Check diagonal deviation constraint
    const std::int64_t diag_dev = std::abs(delta_ref - delta_query);
    if (static_cast<std::size_t>(diag_dev) > config_.max_diag_dev) {
        return false;
    }

    return true;
}

double DPChainClusterer::gap_cost(const Anchor& j, const Anchor& i) const {
    // Compute gap between end of j and start of i
    const std::int64_t j_ref_end = j.ref_coord + static_cast<std::int64_t>(j.length);
    const std::int64_t j_query_end = static_cast<std::int64_t>(j.query_pos + j.length);

    const std::int64_t ref_gap = i.ref_coord - j_ref_end;
    const std::int64_t query_gap = static_cast<std::int64_t>(i.query_pos) - j_query_end;

    // Gap is 0 if there's overlap (negative gap)
    const std::int64_t ref_gap_abs = std::max<std::int64_t>(0, ref_gap);
    const std::int64_t query_gap_abs = std::max<std::int64_t>(0, query_gap);

    // Average gap distance
    const double avg_gap = (ref_gap_abs + query_gap_abs) / 2.0;

    // Distance penalty
    double cost = avg_gap * config_.gap_penalty_factor;

    // Diagonal deviation penalty
    const std::int64_t delta_ref = i.ref_coord - j.ref_coord;
    const std::int64_t delta_query = static_cast<std::int64_t>(i.query_pos) -
                                     static_cast<std::int64_t>(j.query_pos);
    const double diag_dev = std::abs(delta_ref - delta_query);
    cost += diag_dev * config_.diag_penalty_factor;

    // Overlap penalty (if anchors overlap)
    if (ref_gap < 0 || query_gap < 0) {
        const double ref_overlap = std::abs(std::min<std::int64_t>(0, ref_gap));
        const double query_overlap = std::abs(std::min<std::int64_t>(0, query_gap));
        const double avg_overlap = (ref_overlap + query_overlap) / 2.0;
        cost += avg_overlap * config_.overlap_penalty_factor;
    }

    // Path switch penalty
    if (config_.allow_cross_haplotypes && i.path_id != j.path_id) {
        cost += config_.path_switch_cost;
    }

    return cost;
}

double DPChainClusterer::anchor_score(const Anchor& anchor) const {
    // Score based on anchor length (coverage)
    return static_cast<double>(anchor.length) * config_.anchor_weight;
}

std::vector<std::size_t> DPChainClusterer::backtrack_chain(
    const std::vector<int>& pred,
    std::size_t best_idx) const {

    std::vector<std::size_t> chain;

    // Backtrack from best_idx following predecessor pointers
    int current = static_cast<int>(best_idx);
    while (current != -1) {
        chain.push_back(static_cast<std::size_t>(current));
        current = pred[current];
    }

    // Reverse to get forward order (chain starts at beginning of read)
    std::reverse(chain.begin(), chain.end());

    return chain;
}

}  // namespace piru::mapping
