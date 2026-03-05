// SPDX-License-Identifier: MIT
// Factory for chaining backends.

#include "mapping/chainer.hpp"
#include "mapping/dp_chainer.hpp"
#include "mapping/dp_chainer_config.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

ChainerPtr make_chainer(const std::string& backend, const cli::Parsed& parsed) {
    if (backend == "dp-chain") {
        auto config = DPChainerConfig::from_parsed(parsed);
        LOG_DEBUG("DP chain config: max_dist=" + std::to_string(config.max_dist) +
                  ", max_diag_dev=" + std::to_string(config.max_diag_dev) +
                  ", gap_penalty=" + std::to_string(config.gap_penalty_factor) +
                  ", diag_penalty=" + std::to_string(config.diag_penalty_factor) +
                  ", overlap_penalty=" + std::to_string(config.overlap_penalty_factor) +
                  ", anchor_weight=" + std::to_string(config.anchor_weight) +
                  ", min_score=" + std::to_string(config.min_chain_score) +
                  ", max_chains=" + std::to_string(config.max_chains) +
                  ", max_skip=" + std::to_string(config.max_skip) +
                  ", merge_chains=" + (config.merge_chains ? "true" : "false"));
        return std::make_unique<DPChainer>(config);
    }
    throw std::runtime_error("Unknown chainer backend: " + backend);
}

}  // namespace piru::mapping
