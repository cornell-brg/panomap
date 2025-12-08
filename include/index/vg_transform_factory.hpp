// SPDX-License-Identifier: MIT
// Factory for creating VGTransform implementations.

#pragma once

#include "index/vg_transform.hpp"
#include <memory>
#include <string>

namespace piru::index {

// Configuration for VGTransform implementations.
// This struct will be extended as more options for different transform strategies are added.
struct TransformConfig {
    // For PathTraversalTransform: strategy for handling nodes not covered by paths.
    // Options: "ambiguous", "holes", "expand", "skip"
    std::string uncovered_strategy{"expand"};

    // If true, ambiguous contexts (e.g., from averaging k-mer values) will be marked
    // with quality annotations in the output signal.
    bool mark_ambiguous_quality{false};
};

// Factory function to create a unique_ptr to a VGTransform implementation.
// backend: Specifies which transformation strategy to use (e.g., "expansion", "path_traversal").
// config: Configuration specific to the chosen transformation strategy.
std::unique_ptr<VGTransform> makeVGTransform(const std::string& backend, const TransformConfig& config);

} // namespace piru::index
