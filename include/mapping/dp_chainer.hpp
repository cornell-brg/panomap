// SPDX-License-Identifier: MIT
// DP-based colinear chainer for path-aware seed selection.

#pragma once

#include <cstddef>
#include <vector>

#include "cli/parse.hpp"
#include "index/linearizer.hpp"
#include "mapping/anchor_expander.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

// Configuration for DP-based chaining algorithm.
// Defaults tuned for noisy nanopore signals (DEV027).
struct DPChainerConfig {
    std::size_t max_dist{500};        // Max query/ref distance for chaining (banding)
    std::size_t max_diag_dev{500};    // Max diagonal deviation |dr - dq|
    std::size_t min_chain_score{12};  // Min score to report a chain
    std::size_t max_chains{10};       // Max number of chains to extract (multi-mapping)
    std::size_t max_skip{25};         // Stop after N consecutive failed chain attempts

    double anchor_weight{1.0};            // Weight per anchor length
    double gap_penalty_factor{0.02};      // Penalty per unit gap distance
    double diag_penalty_factor{0.05};     // Penalty per unit diagonal deviation
    double overlap_penalty_factor{0.90};  // Penalty per unit overlap

    bool merge_chains{false};  // Merge overlapping chains on same path

    // CLI integration: options and parsing for --chain-* flags.
    static std::vector<cli::Option> cli_options();
    static DPChainerConfig from_parsed(const cli::Parsed& parsed);
};

// DP-based colinear chainer.
//
// Algorithm:
// 1. Sort anchors by (path_id, ref_coord, query_pos)
// 2. DP: dp[i] = max_j(dp[j] + score(i) - gap_cost(j, i))
//    where j is a valid predecessor (same path, within distance, diagonal constraints)
// 3. Backtrack from best scoring anchor to extract chain
// 4. Optionally extract multiple chains for multi-mapping
//
// Returns ChainResult with anchors from extracted chain(s).
//
// Note: Operates on anchors (linear space). Expansion from seed hits must be done
// by caller using AnchorExpander.
class DPChainer : public Chainer {
public:
    // Construct with config (no longer needs linearization_coords).
    explicit DPChainer(DPChainerConfig config);

    ChainResult chain(const std::vector<PathAnchor>& anchors) const override;
    std::string name() const override { return "dp-chain"; }

    void dump_path_chains(const char* filename, const std::string& read_id,
                          std::size_t read_length, const std::vector<PathAnchor>& anchors,
                          const index::GraphStore* graph_store) const override;

    void dump_anchor_detail(const char* filename, const std::string& read_id,
                            std::size_t read_length, const std::vector<PathAnchor>& anchors,
                            const index::GraphStore* graph_store) const override;

private:
    DPChainerConfig config_;

    // Check if anchor j can chain to anchor i.
    bool can_chain(const PathAnchor& j, const PathAnchor& i) const;

    // Compute gap cost for chaining anchor j to anchor i.
    double gap_cost(const PathAnchor& j, const PathAnchor& i) const;

    // Compute anchor score (based on length).
    double anchor_score(const PathAnchor& anchor) const;

    // Extract chain by backtracking from given anchor index.
    std::vector<std::size_t> backtrack_chain(const std::vector<int>& pred,
                                             std::size_t best_idx) const;
};

}  // namespace piru::mapping
