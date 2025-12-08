// SPDX-License-Identifier: MIT
// Functions for serializing and deserializing index components.

#pragma once

#include <string>
#include <memory>
#include <utility>

#include "index/graph_store.hpp"
#include "index/signal_store.hpp"

namespace piru::io::index {

// Corresponds to the Global Index Metadata section in the .graph file
struct IndexMetadata {
    uint32_t piru_version_major{0};
    uint32_t piru_version_minor{0};
    uint32_t piru_version_patch{0};
    uint64_t build_timestamp{0};
    uint32_t graph_flavor{0}; // 0=unknown, 1=dbg, 2=vg
    uint32_t graph_k{0};
    uint32_t pore_k{0};
    std::string model_name;
    std::string fuzzy_quantizer;
    std::string align_quantizer;
    std::string source_path;
};

// Writes the GraphStore to the specified file path.
void write_graph(const std::string& path,
                 const piru::index::GraphStore& store,
                 const IndexMetadata& metadata);

// Reads a GraphStore from the specified file path.
// Returns the loaded graph store and the global metadata read from the file.
std::pair<std::unique_ptr<piru::index::AdjListGraphStore>, IndexMetadata>
read_graph(const std::string& path);

struct SignalMetadata {
    uint32_t quantization_bits{0};
    float scale{1.0f};
    float offset{0.0f};
};

// Writes the SignalStore to the specified file path.
void write_signals(const std::string& path, const piru::index::SignalStore& store, float scale, float offset);

// Reads a SignalStore from the specified file path.
std::pair<std::unique_ptr<piru::index::VectorSignalStore>, SignalMetadata>
read_signals(const std::string& path);

} // namespace piru::io::index
