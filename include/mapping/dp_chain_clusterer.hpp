// SPDX-License-Identifier: MIT
// DP-based colinear chaining clusterer for path-aware seed selection.

#pragma once

#include <cstddef>
#include <vector>

#include "index/linearizer.hpp"
#include "mapping/anchor_expander.hpp"
#include "mapping/seed_clusterer.hpp"

namespace piru::mapping {

// Configuration for DP-based chaining algorithm.
struct DPChainClustererConfig {
    std::size_t max_dist{5000};              // Max query/ref distance for chaining
    std::size_t max_diag_dev{500};           // Max diagonal deviation |Δr - Δq|
    bool allow_cross_haplotypes{false};      // Allow chaining across different paths
    double path_switch_cost{50.0};           // Penalty for switching between paths
    std::size_t min_chain_score{100};        // Min score to report a chain
    std::size_t max_chains{10};              // Max number of chains to extract (multi-mapping)

    // Scoring parameters
    double anchor_weight{1.0};               // Weight per anchor length
    double gap_penalty_factor{0.1};          // Penalty per unit gap distance
    double diag_penalty_factor{0.5};         // Penalty per unit diagonal deviation
    double overlap_penalty_factor{2.0};      // Penalty per unit overlap
};

// DP-based colinear chaining clusterer.
//
// Algorithm:
// 1. Expand seed hits to anchors using linearization coordinates
// 2. Sort anchors by (path_id, ref_coord, query_pos)
// 3. DP: dp[i] = max_j(dp[j] + score(i) - gap_cost(j, i))
//    where j is a valid predecessor (same/allowed path, within distance, diagonal constraints)
// 4. Backtrack from best scoring anchor to extract chain
// 5. Optionally extract multiple chains for multi-mapping
//
// Returns ClusterSummary with anchors from extracted chain(s).
class DPChainClusterer : public SeedClusterer {
public:
    // Construct with linearization coordinates and config.
    DPChainClusterer(const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     DPChainClustererConfig config);

    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override;
    std::string name() const override { return "dp-chain"; }

private:
    AnchorExpander expander_;
    DPChainClustererConfig config_;

    // Check if anchor j can chain to anchor i.
    bool can_chain(const Anchor& j, const Anchor& i) const;

    // Compute gap cost for chaining anchor j to anchor i.
    double gap_cost(const Anchor& j, const Anchor& i) const;

    // Compute anchor score (based on length).
    double anchor_score(const Anchor& anchor) const;

    // Extract chain by backtracking from given anchor index.
    std::vector<std::size_t> backtrack_chain(const std::vector<int>& pred,
                                              std::size_t best_idx) const;
};

}  // namespace piru::mapping
