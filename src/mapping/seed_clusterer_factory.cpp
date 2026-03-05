// SPDX-License-Identifier: MIT
// Factory for seed clustering/chaining backends.

#include <stdexcept>

#include "mapping/dp_chain_clusterer.hpp"
#include "mapping/seed_clusterer.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

AnchorClustererPtr make_anchor_clusterer(const SeedClustererConfig& config) {
    if (config.backend == "dp-chain" || config.backend.empty()) {
        DPChainClustererConfig dp_config;
        dp_config.max_dist = config.dp_max_dist;
        dp_config.max_diag_dev = config.dp_max_diag_dev;
        dp_config.gap_penalty_factor = config.dp_gap_penalty;
        dp_config.diag_penalty_factor = config.dp_diag_penalty;
        dp_config.overlap_penalty_factor = config.dp_overlap_penalty;
        dp_config.anchor_weight = config.dp_anchor_weight;
        dp_config.min_chain_score = config.dp_min_chain_score;
        dp_config.max_chains = config.dp_max_chains;
        dp_config.max_skip = config.dp_max_skip;
        dp_config.merge_chains = config.dp_merge_chains;
        LOG_DEBUG("DP chain config: max_dist=" + std::to_string(dp_config.max_dist) +
                  ", max_diag_dev=" + std::to_string(dp_config.max_diag_dev) +
                  ", gap_penalty=" + std::to_string(dp_config.gap_penalty_factor) +
                  ", diag_penalty=" + std::to_string(dp_config.diag_penalty_factor) +
                  ", overlap_penalty=" + std::to_string(dp_config.overlap_penalty_factor) +
                  ", anchor_weight=" + std::to_string(dp_config.anchor_weight) +
                  ", min_score=" + std::to_string(dp_config.min_chain_score) +
                  ", max_chains=" + std::to_string(dp_config.max_chains) +
                  ", max_skip=" + std::to_string(dp_config.max_skip) +
                  ", merge_chains=" + (dp_config.merge_chains ? "true" : "false"));
        return std::make_unique<DPChainClusterer>(dp_config);
    }
    throw std::runtime_error("Unknown clusterer backend: " + config.backend);
}

}  // namespace piru::mapping
