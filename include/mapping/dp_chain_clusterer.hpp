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
// Note: Cross-path chaining is not supported. See DEV033 for planned alias-based approach.
struct DPChainClustererConfig {
    std::size_t max_dist{500};        // Max query/ref distance for chaining (also used for banding)
    std::size_t max_diag_dev{500};    // Max diagonal deviation |Δr - Δq|
    std::size_t min_chain_score{12};  // Min score to report a chain
    std::size_t max_chains{10};       // Max number of chains to extract (multi-mapping)
    std::size_t max_skip{25};         // Stop after this many consecutive failed chain attempts

    // Scoring parameters
    double anchor_weight{1.5};           // Weight per anchor length
    double gap_penalty_factor{0.1};      // Penalty per unit gap distance
    double diag_penalty_factor{0.5};     // Penalty per unit diagonal deviation
    double overlap_penalty_factor{0.4};  // Penalty per unit overlap

    // Post-processing
    bool merge_chains{true};  // Merge overlapping chains on same path
};

// DP-based colinear chaining clusterer.
//
// Algorithm:
// 1. Sort anchors by (path_id, ref_coord, query_pos)
// 2. DP: dp[i] = max_j(dp[j] + score(i) - gap_cost(j, i))
//    where j is a valid predecessor (same path, within distance, diagonal constraints)
// 3. Backtrack from best scoring anchor to extract chain
// 4. Optionally extract multiple chains for multi-mapping
//
// Returns ClusterSummary with anchors from extracted chain(s).
//
// Note: Operates on anchors (linear space). Expansion from seed hits must be done
// by caller using AnchorExpander.
class DPChainClusterer : public AnchorClusterer {
public:
    // Construct with config (no longer needs linearization_coords).
    explicit DPChainClusterer(DPChainClustererConfig config);

    ClusterSummary cluster(const std::vector<Anchor>& anchors) const override;
    std::string name() const override { return "dp-chain"; }

private:
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
