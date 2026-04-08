/**
 * serialization.hpp
 *
 * .pirx index format: single-file binary containing graph topology,
 * linearization coordinates, and seed hash table.
 *
 * Format: PIRX magic + metadata + graph (nodes/edges/paths)
 *         + linearization + seeds.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"

namespace piru::io::index {

struct IndexMetadata {
  std::string version;            // piru version that built this index
  uint64_t build_timestamp{0};    // unix seconds
  std::string model_name;
  uint32_t pore_k{0};
  std::string tokenizer;
};

struct LoadedIndex {
  IndexMetadata metadata;
  std::unique_ptr<piru::index::FlatGraphStore> graph;
  std::unique_ptr<piru::index::SeedStore> seeds;
  std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords;
  std::vector<float> node_1d_coords;          // empty if not computed at index time
  std::vector<std::uint32_t> component_ids;  // connected component per node
};

/** Save an index to a single .pirx file. */
void save_index(const std::string& path, const piru::index::GraphStore& graph,
                const piru::index::SeedStore& seeds,
                const std::vector<std::vector<piru::index::LinearCoordinate>>& linearization_coords,
                const IndexMetadata& metadata, const std::vector<float>& node_1d_coords = {},
                const std::vector<std::uint32_t>& component_ids = {});

/** Load an index from a .pirx file. */
LoadedIndex load_index(const std::string& path);

/** Check if a file has the PIRX magic header. */
bool is_pirx_index(const std::string& path);

}  // namespace piru::io::index
