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

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "core/index/graph_store.hpp"
#include "core/index/linearizer.hpp"
#include "core/index/seed_store.hpp"

namespace panomap::io::index {

// Index modality. Stamped in the .pirx header by the binary that built it,
// checked at load time so a panomap binary refuses a panomap-base index
// (and vice versa). Format break: v1 indexes have no mode byte and are
// rejected by v2+ loaders.
enum class IndexMode : std::uint8_t {
  kSignal = 0,
  kBase = 1,
};

const char* mode_name(IndexMode mode);

struct IndexMetadata {
  std::string version;            // panomap version that built this index
  uint64_t build_timestamp{0};    // unix seconds
  IndexMode mode{IndexMode::kSignal};
  std::string model_name;
  uint32_t pore_k{0};
  std::string tokenizer;
};

// Byte sizes of each section in the .pirx file. Populated by load_index;
// used by `panomap inspect` to show on-disk breakdown.
struct SectionSizes {
  uint64_t header_meta{0};   // magic + version + flags + mode + metadata strings
  uint64_t graph{0};         // nodes + edges + paths
  uint64_t linearization{0}; // per-node (path_id, ref_coord) array
  uint64_t seeds{0};         // bucket-native seed store (keys + counts + offsets + entries)
  uint64_t coords_1d{0};     // float32 per node (optional)
  uint64_t components{0};    // uint32 per node (optional)
  uint64_t total{0};
};

struct LoadedIndex {
  IndexMetadata metadata;
  std::unique_ptr<panomap::index::FlatGraphStore> graph;
  std::unique_ptr<panomap::index::SeedStore> seeds;
  std::vector<std::vector<panomap::index::LinearCoordinate>> linearization_coords;
  std::vector<float> node_1d_coords;          // empty if not computed at index time
  std::vector<std::uint32_t> component_ids;  // connected component per node
  SectionSizes section_sizes;
};

/** Save an index to a single .pirx file. */
void save_index(const std::string& path, const panomap::index::GraphStore& graph,
                const panomap::index::SeedStore& seeds,
                const std::vector<std::vector<panomap::index::LinearCoordinate>>& linearization_coords,
                const IndexMetadata& metadata, const std::vector<float>& node_1d_coords = {},
                const std::vector<std::uint32_t>& component_ids = {});

/** Load an index from a .pirx file. */
LoadedIndex load_index(const std::string& path);

/** Check if a file has the PIRX magic header. */
bool is_pirx_index(const std::string& path);

}  // namespace panomap::io::index
