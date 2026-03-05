// SPDX-License-Identifier: MIT
// Serialization for .pirx index format (single-file binary).

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"

namespace piru::io::index {

struct IndexMetadata {
    std::string model_name;
    uint32_t pore_k{0};
    std::string fuzzy_quantizer;
};

struct LoadedIndex {
    IndexMetadata metadata;
    std::unique_ptr<piru::index::AdjListGraphStore> graph;
    std::unique_ptr<piru::index::HashSeedStore> seeds;
    std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
};

/**
 * Save an index to a single .pirx file.
 *
 * Format: PIRX header + metadata + graph + linearization + seeds
 */
void save_index(const std::string& path, const piru::index::GraphStore& graph,
                const piru::index::SeedStore& seeds,
                const std::vector<std::vector<piru::index::LinearCoordinate>>& linearization_coords,
                const IndexMetadata& metadata);

/**
 * Load an index from a .pirx file.
 */
LoadedIndex load_index(const std::string& path);

/**
 * Check if a file is a .pirx index (has PIRX magic header).
 */
bool is_pirx_index(const std::string& path);

}  // namespace piru::io::index
