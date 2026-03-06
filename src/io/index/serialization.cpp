/**
 * serialization.cpp
 *
 * Binary read/write for the .pirx index format.
 * See serialization.hpp for format layout.
 *
 * SPDX-License-Identifier: MIT
 */

#include "io/index/serialization.hpp"

#include <fstream>
#include <stdexcept>

#include "index/aln_graph.hpp"
#include "util/logging.hpp"

namespace piru::io::index {

namespace {

constexpr char kMagic[4] = {'P', 'I', 'R', 'X'};
constexpr uint32_t kVersion = 1;

constexpr uint32_t kFlagHasSequences = 1 << 0;
constexpr uint32_t kFlagHasFuzzySignals = 1 << 1;   // reserved
constexpr uint32_t kFlagHasAlignSignals = 1 << 2;    // reserved

template <typename T>
void write_pod(std::ostream& out, const T& val) {
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
void read_pod(std::istream& in, T& val) {
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
}

void write_string(std::ostream& out, const std::string& s) {
    uint32_t len = static_cast<uint32_t>(s.size());
    write_pod(out, len);
    out.write(s.data(), len);
}

std::string read_string(std::istream& in) {
    uint32_t len = 0;
    read_pod(in, len);
    std::string s(len, '\0');
    in.read(s.data(), len);
    return s;
}

}  // namespace

void save_index(const std::string& path, const piru::index::GraphStore& graph_store,
                const piru::index::SeedStore& seed_store,
                const std::vector<std::vector<piru::index::LinearCoordinate>>& linearization_coords,
                const IndexMetadata& metadata) {
    const auto* adj_store = dynamic_cast<const piru::index::AdjListGraphStore*>(&graph_store);
    if (!adj_store) {
        throw std::runtime_error("Unsupported GraphStore backend for serialization");
    }

    const auto* hash_store = dynamic_cast<const piru::index::HashSeedStore*>(&seed_store);
    if (!hash_store) {
        throw std::runtime_error("Unsupported SeedStore backend for serialization");
    }

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    const auto& graph = adj_store->graph();

    /* 1. Header */

    out.write(kMagic, 4);
    write_pod<uint32_t>(out, kVersion);
    uint32_t flags = 0;
    write_pod<uint32_t>(out, flags);

    /* 2. Metadata */

    write_string(out, metadata.model_name);
    write_pod<uint32_t>(out, metadata.pore_k);
    write_string(out, metadata.fuzzy_quantizer);
    // Legacy: graph_flavor was "vg" vs "dbg", now unused. Write empty string
    // to preserve .pirx binary layout. Remove on next format version bump.
    write_string(out, std::string{});

    /* 3. Graph - nodes */

    write_pod<uint64_t>(out, graph.nodeCount());
    for (std::size_t i = 0; i < graph.nodeCount(); ++i) {
        const auto& node = graph.node(i);
        write_string(out, node.original_id);
        write_pod<uint8_t>(out, node.is_reverse ? 1 : 0);
        if (flags & kFlagHasSequences) {
            write_string(out, node.sequence);
        }
    }

    /* 4. Graph - edges */

    write_pod<uint64_t>(out, graph.edgeCount());
    for (const auto& edge : graph.edges()) {
        write_pod<uint64_t>(out, edge.from);
        write_pod<uint64_t>(out, edge.to);
        write_pod<uint64_t>(out, edge.overlap_bases);
    }

    /* 5. Graph - paths */

    write_pod<uint64_t>(out, graph.pathCount());
    for (std::size_t i = 0; i < graph.pathCount(); ++i) {
        const auto& path = graph.paths()[i];
        write_string(out, path.name);
        write_pod<uint64_t>(out, path.length);
        write_pod<uint64_t>(out, path.steps.size());
        for (const auto& step : path.steps) {
            uint64_t node_id = std::stoull(step.node_id);
            write_pod<uint64_t>(out, node_id);
        }
    }

    /* 6. Linearization */

    write_pod<uint64_t>(out, linearization_coords.size());
    for (const auto& coords : linearization_coords) {
        write_pod<uint64_t>(out, coords.size());
        for (const auto& coord : coords) {
            write_pod<uint64_t>(out, coord.path_id);
            write_pod<int64_t>(out, coord.ref_coord);
        }
    }

    /* 7. Seeds */

    write_string(out, hash_store->extractor_name());
    const auto& params = hash_store->params();
    write_pod<uint64_t>(out, params.size());
    for (const auto& [key, value] : params) {
        write_string(out, key);
        write_string(out, value);
    }
    write_pod<double>(out, hash_store->filter_fraction());
    write_pod<uint64_t>(out, hash_store->max_hash_frequency());
    write_pod<uint64_t>(out, hash_store->frequency_threshold());

    const auto& seed_data = hash_store->data();
    write_pod<uint64_t>(out, seed_data.size());
    for (const auto& [hash, hits] : seed_data) {
        write_pod<uint64_t>(out, hash);
        write_pod<uint64_t>(out, hits.size());
        for (const auto& hit : hits) {
            write_pod<uint64_t>(out, hit.node_id);
            write_pod<uint64_t>(out, hit.offset);
            write_pod<uint64_t>(out, hit.length);
        }
    }

    LOG_INFO("Saved index: " + std::to_string(graph.nodeCount()) + " nodes, " +
             std::to_string(seed_data.size()) + " seeds -> " + path);
}

LoadedIndex load_index(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    /* 1. Header */

    char magic[4];
    in.read(magic, 4);
    if (std::string(magic, 4) != std::string(kMagic, 4)) {
        throw std::runtime_error("Invalid magic number for index file");
    }

    uint32_t version;
    read_pod(in, version);
    if (version > kVersion) {
        LOG_WARN("Index version " + std::to_string(version) + " is newer than supported version " +
                 std::to_string(kVersion));
    }

    uint32_t flags;
    read_pod(in, flags);

    /* 2. Metadata */

    IndexMetadata metadata;
    metadata.model_name = read_string(in);
    read_pod(in, metadata.pore_k);
    metadata.fuzzy_quantizer = read_string(in);
    // Legacy: graph_flavor, read and discard
    read_string(in);

    /* 3. Graph - nodes */

    uint64_t node_count;
    read_pod(in, node_count);

    auto graph = std::make_unique<piru::index::AlnGraph>();
    graph->reserveNodes(node_count);

    for (uint64_t i = 0; i < node_count; ++i) {
        piru::index::AlnNode node;
        node.original_id = read_string(in);
        uint8_t is_reverse;
        read_pod(in, is_reverse);
        node.is_reverse = (is_reverse != 0);
        if (flags & kFlagHasSequences) {
            node.sequence = read_string(in);
        }
        graph->setNode(i, std::move(node));
    }

    /* 4. Graph - edges */

    uint64_t edge_count;
    read_pod(in, edge_count);
    for (uint64_t i = 0; i < edge_count; ++i) {
        piru::index::AlnEdge edge;
        read_pod(in, edge.from);
        read_pod(in, edge.to);
        read_pod(in, edge.overlap_bases);
        graph->addEdge(edge);
    }

    /* 5. Graph - paths */

    uint64_t path_count;
    read_pod(in, path_count);
    for (uint64_t i = 0; i < path_count; ++i) {
        piru::index::AlnPath aln_path;
        aln_path.name = read_string(in);
        read_pod(in, aln_path.length);

        uint64_t step_count;
        read_pod(in, step_count);
        aln_path.steps.reserve(step_count);
        for (uint64_t j = 0; j < step_count; ++j) {
            uint64_t node_id;
            read_pod(in, node_id);
            piru::index::AlnPathStep step;
            step.node_id = std::to_string(node_id);
            aln_path.steps.push_back(std::move(step));
        }
        graph->addPath(std::move(aln_path));
    }

    /* 6. Linearization */

    uint64_t lin_node_count;
    read_pod(in, lin_node_count);

    std::vector<std::vector<piru::index::LinearCoordinate>> linearization_coords(lin_node_count);
    for (uint64_t i = 0; i < lin_node_count; ++i) {
        uint64_t coord_count;
        read_pod(in, coord_count);
        linearization_coords[i].reserve(coord_count);
        for (uint64_t j = 0; j < coord_count; ++j) {
            uint64_t path_id;
            int64_t ref_coord;
            read_pod(in, path_id);
            read_pod(in, ref_coord);
            linearization_coords[i].emplace_back(path_id, ref_coord);
        }
    }

    /* 7. Seeds */

    auto seeds = std::make_unique<piru::index::HashSeedStore>();
    seeds->set_extractor_name(read_string(in));

    uint64_t param_count;
    read_pod(in, param_count);
    std::map<std::string, std::string> params;
    for (uint64_t i = 0; i < param_count; ++i) {
        std::string key = read_string(in);
        std::string value = read_string(in);
        params[key] = value;
    }
    seeds->set_params(std::move(params));

    double filter_fraction;
    read_pod(in, filter_fraction);
    seeds->set_filter_fraction(filter_fraction);

    uint64_t max_freq;
    read_pod(in, max_freq);
    seeds->set_max_hash_frequency(max_freq);

    uint64_t freq_threshold;
    read_pod(in, freq_threshold);
    seeds->set_frequency_threshold(freq_threshold);

    uint64_t entry_count;
    read_pod(in, entry_count);
    for (uint64_t i = 0; i < entry_count; ++i) {
        uint64_t hash;
        read_pod(in, hash);

        uint64_t hit_count;
        read_pod(in, hit_count);
        for (uint64_t j = 0; j < hit_count; ++j) {
            uint64_t node_id, offset, length;
            read_pod(in, node_id);
            read_pod(in, offset);
            read_pod(in, length);
            seeds->insert(hash, piru::index::SeedEntry{node_id, offset, length});
        }
    }

    LOG_INFO("Loaded index: " + std::to_string(node_count) + " nodes, " +
             std::to_string(entry_count) + " seeds ← " + path);

    return {std::move(metadata),
            std::make_unique<piru::index::AdjListGraphStore>(std::move(*graph)), std::move(seeds),
            std::move(linearization_coords)};
}

bool is_pirx_index(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    char magic[4];
    in.read(magic, 4);
    return std::string(magic, 4) == std::string(kMagic, 4);
}

}  // namespace piru::io::index
