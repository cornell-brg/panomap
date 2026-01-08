// SPDX-License-Identifier: MIT
// Serialization for simple index pipeline (single-file format).

#pragma once

#include <string>
#include <memory>
#include <vector>

#include "index/graph_store.hpp"
#include "index/seed_store.hpp"
#include "index/linearizer.hpp"

namespace piru::io::index {

struct SimpleIndexMetadata {
    std::string model_name;
    uint32_t pore_k{0};
    std::string fuzzy_quantizer;
    std::string graph_flavor;
};

struct SimpleLoadedIndex {
    SimpleIndexMetadata metadata;
    std::unique_ptr<piru::index::AdjListGraphStore> graph;
    std::unique_ptr<piru::index::HashSeedStore> seeds;
    std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
};

/**
 * Save a simple index to a single file.
 *
 * Format: PIR2 header + metadata + graph + linearization + seeds
 */
void save_simple_index(
    const std::string& path,
    const piru::index::GraphStore& graph,
    const piru::index::SeedStore& seeds,
    const std::vector<std::vector<piru::index::LinearCoordinate>>& linearization_coords,
    const SimpleIndexMetadata& metadata);

/**
 * Load a simple index from a single file.
 */
SimpleLoadedIndex load_simple_index(const std::string& path);

/**
 * Check if a file is a simple index (has PIR2 magic).
 */
bool is_simple_index(const std::string& path);

}  // namespace piru::io::index
