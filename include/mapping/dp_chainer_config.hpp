// SPDX-License-Identifier: MIT
// Configuration for the DP-based colinear chainer.

#pragma once

#include <cstddef>
#include <vector>

#include "cli/parse.hpp"

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

}  // namespace piru::mapping
