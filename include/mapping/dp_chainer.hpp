// SPDX-License-Identifier: MIT
// DP-based colinear chainer for path-aware seed selection.

#pragma once

#include <cstddef>
#include <vector>

#include "index/linearizer.hpp"
#include "mapping/anchor_expander.hpp"
#include "mapping/dp_chainer_config.hpp"
#include "mapping/chainer.hpp"

namespace piru::mapping {

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

    ChainResult chain(const std::vector<Anchor>& anchors) const override;
    std::string name() const override { return "dp-chain"; }

    void dump_path_chains(const char* filename, const std::string& read_id,
                          std::size_t read_length, const std::vector<Anchor>& anchors,
                          const index::GraphStore* graph_store) const override;

    void dump_anchor_detail(const char* filename, const std::string& read_id,
                            std::size_t read_length, const std::vector<Anchor>& anchors,
                            const index::GraphStore* graph_store) const override;

private:
    DPChainerConfig config_;

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
