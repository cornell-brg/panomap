// SPDX-License-Identifier: MIT
#include "io/index/serialization.hpp"

#include <chrono>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "version.hpp"
#include "index/aln_graph.hpp"
#include "index/signal_store.hpp"
#include "index/seed_store.hpp"
#include "util/logging.hpp"

namespace piru::io::index {
namespace {

// Helper to write a plain old data (POD) type to a stream.
template <typename T>
void write_pod(std::ostream& out, const T& val) {
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

// Helper to read a plain old data (POD) type from a stream.
template <typename T>
void read_pod(std::istream& in, T& val) {
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
}

// Helper to write a std::string to a stream (length-prefixed).
void write_string(std::ostream& out, const std::string& s) {
    uint32_t len = static_cast<uint32_t>(s.length());
    write_pod(out, len);
    out.write(s.c_str(), len);
}

// Helper to read a std::string from a stream (length-prefixed).
std::string read_string(std::istream& in) {
    uint32_t len = 0;
    read_pod(in, len);
    std::vector<char> buf(len);
    in.read(buf.data(), len);
    return std::string(buf.begin(), buf.end());
}

} // namespace

void write_graph(const std::string& path,
                 const piru::index::GraphStore& store,
                 const IndexMetadata& metadata) {
    // --- Backend detection ---
    const auto* adj_store = dynamic_cast<const piru::index::AdjListGraphStore*>(&store);
    if (!adj_store) {
        throw std::runtime_error("Unsupported GraphStore backend type for serialization.");
    }
    const std::string backend_type = "adjlist";

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    // --- Common Header ---
    const char magic[8] = {'P', 'I', 'R', 'U', 'G', 'R', 'A', 'F'};
    out.write(magic, 8);
    write_pod<uint32_t>(out, 1000); // Format version 1.0

    // Placeholder for header size, we'll come back and write it later.
    const auto header_size_pos = out.tellp();
    write_pod<uint32_t>(out, 0);

    write_string(out, backend_type);

    const auto metadata_start_pos = out.tellp();

    // --- Global Index Metadata ---
    write_pod<uint32_t>(out, piru::Version::kMajor);
    write_pod<uint32_t>(out, piru::Version::kMinor);
    write_pod<uint32_t>(out, piru::Version::kPatch);

    const auto timestamp = std::chrono::system_clock::now();
    const auto timestamp_since_epoch = std::chrono::duration_cast<std::chrono::seconds>(
        timestamp.time_since_epoch()).count();
    write_pod<uint64_t>(out, timestamp_since_epoch);

    write_pod<uint32_t>(out, metadata.graph_flavor);
    write_pod<uint32_t>(out, metadata.graph_k);
    write_pod<uint32_t>(out, metadata.pore_k);
    write_string(out, metadata.model_name);
    write_string(out, metadata.fuzzy_quantizer);
    write_string(out, metadata.align_quantizer);
    write_string(out, metadata.source_path);

    // --- Graph Metadata ---
    const auto& graph = adj_store->graph();
    write_pod<uint64_t>(out, graph.nodeCount());
    write_pod<uint64_t>(out, graph.edgeCount());
    write_pod<uint64_t>(out, 0); // Path count, not supported yet
    write_pod<uint32_t>(out, 0); // Reserved
    write_pod<uint32_t>(out, 0); // Reserved

    // Now calculate and write the actual header size.
    const auto data_start_pos = out.tellp();
    const uint32_t header_size = static_cast<uint32_t>(data_start_pos);
    out.seekp(header_size_pos);
    write_pod<uint32_t>(out, header_size);
    out.seekp(data_start_pos);

    // --- Node Array ---
    for (size_t i = 0; i < graph.nodeCount(); ++i) {
        const auto& node = graph.node(i);
        write_string(out, node.sequence);
        write_pod<uint32_t>(out, node.chain_id.value_or(-1));
        write_pod<int64_t>(out, node.linear_position.value_or(0));
        write_string(out, node.label);
    }

    // --- Edge Array ---
    for (const auto& edge : graph.edges()) {
        write_pod<uint32_t>(out, edge.from);
        write_pod<uint32_t>(out, edge.to);
        write_pod<uint32_t>(out, edge.overlap_bases);
    }

    // --- Path Array ---
    // Not implemented.
}

std::pair<std::unique_ptr<piru::index::AdjListGraphStore>, IndexMetadata>
read_graph(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    // --- Common Header ---
    char magic[8];
    in.read(magic, 8);
    if (std::string(magic, 8) != "PIRUGRAF") {
        throw std::runtime_error("Invalid magic number for .graph file.");
    }

    uint32_t format_version;
    read_pod(in, format_version);
    if (format_version > 1000) {
        LOG_WARN("Attempting to read a .graph file with a newer version (" +
                 std::to_string(format_version) + "). Compatibility is not guaranteed.");
    }

    uint32_t header_size;
    read_pod(in, header_size);

    std::string backend_type = read_string(in);
    if (backend_type != "adjlist") {
        throw std::runtime_error("Unsupported backend type for .graph file: " + backend_type);
    }

    // --- Global Index Metadata ---
    IndexMetadata metadata;
    read_pod(in, metadata.piru_version_major);
    read_pod(in, metadata.piru_version_minor);
    read_pod(in, metadata.piru_version_patch);
    read_pod(in, metadata.build_timestamp);
    read_pod(in, metadata.graph_flavor);
    read_pod(in, metadata.graph_k);
    read_pod(in, metadata.pore_k);
    metadata.model_name = read_string(in);
    metadata.fuzzy_quantizer = read_string(in);
    metadata.align_quantizer = read_string(in);
    metadata.source_path = read_string(in);

    // --- Graph Metadata ---
    uint64_t node_count, edge_count, path_count;
    read_pod(in, node_count);
    read_pod(in, edge_count);
    read_pod(in, path_count);
    in.seekg(8, std::ios_base::cur); // Skip reserved fields

    if (static_cast<uint32_t>(in.tellg()) != header_size) {
        throw std::runtime_error("Header size mismatch in .graph file.");
    }

    auto graph = std::make_unique<piru::index::AlnGraph>();
    // --- Node Array ---
    for (uint64_t i = 0; i < node_count; ++i) {
        piru::index::AlnNode node;
        node.sequence = read_string(in);
        uint32_t chain_id_u32;
        read_pod(in, chain_id_u32);
        if (chain_id_u32 != static_cast<uint32_t>(-1)) {
            node.chain_id = chain_id_u32;
        }
        int64_t linear_pos;
        read_pod(in, linear_pos);
        node.linear_position = linear_pos;
        node.label = read_string(in);
        graph->addNode(std::move(node));
    }

    // --- Edge Array ---
    for (uint64_t i = 0; i < edge_count; ++i) {
        piru::index::AlnEdge edge;
        uint32_t from, to, overlap;
        read_pod(in, from);
        read_pod(in, to);
        read_pod(in, overlap);
        edge.from = from;
        edge.to = to;
        edge.overlap_bases = overlap;
        graph->addEdge(edge);
    }
    
    // --- Path Array ---
    // Not implemented in writer, but need to handle skipping if present.
    if(path_count > 0) {
        LOG_WARN("Path array found in .graph file but is not supported. Skipping.");
        // This is a placeholder for skipping logic. For now, we assume it's not there.
    }

    // Verify that we have reached the end of the file.
    if (in.peek() != EOF) {
        throw std::runtime_error("Trailing data found in .graph file after expected content.");
    }

    auto store = std::make_unique<piru::index::AdjListGraphStore>(std::move(*graph));

    return {std::move(store), std::move(metadata)};
}


void write_signals(const std::string& path, const piru::index::SignalStore& store, float scale, float offset) {
    // --- Backend detection ---
    const auto* vec_store = dynamic_cast<const piru::index::VectorSignalStore*>(&store);
    if (!vec_store) {
        throw std::runtime_error("Unsupported SignalStore backend type for serialization.");
    }

    const auto& signals = vec_store->signals();
    std::string backend_type = "float32"; // Default
    uint32_t quantization_bits = 32;

    if (!signals.empty()) {
        const auto first_kind = signals[0].kind;
        if (first_kind == piru::signal::AlignmentQuantizationKind::kInt16) {
            backend_type = "int16";
            quantization_bits = 16;
        } else if (first_kind == piru::signal::AlignmentQuantizationKind::kInt8) {
            backend_type = "int8";
            quantization_bits = 8;
        } else if (first_kind != piru::signal::AlignmentQuantizationKind::kFloat32) {
             throw std::runtime_error("Unsupported signal quantization kind for serialization.");
        }
        
        // Validate that all signals have the same kind.
        for(size_t i = 1; i < signals.size(); ++i) {
            if (signals[i].kind != first_kind) {
                throw std::runtime_error("Inconsistent signal kinds in SignalStore.");
            }
        }
    }
    
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    // --- Common Header ---
    const char magic[8] = {'P', 'I', 'R', 'U', 'S', 'I', 'G', 'S'};
    out.write(magic, 8);
    write_pod<uint32_t>(out, 1000); // Format version 1.0

    const auto header_size_pos = out.tellp();
    write_pod<uint32_t>(out, 0); // Placeholder for header size
    
    write_string(out, backend_type);

    // --- Signal Metadata ---
    write_pod<uint64_t>(out, signals.size());
    write_pod<uint32_t>(out, quantization_bits);
    write_pod<float>(out, scale);
    write_pod<float>(out, offset);
    write_pod<uint32_t>(out, 0);  // Reserved

    const uint32_t header_size = static_cast<uint32_t>(out.tellp());
    out.seekp(header_size_pos);
    write_pod<uint32_t>(out, header_size);
    out.seekp(header_size);

    // --- Signal Array ---
    for (const auto& signal : signals) {
        std::visit([&out](const auto& data_vec) {
            using T = typename std::decay_t<decltype(data_vec)>::value_type;
            uint32_t sample_count = static_cast<uint32_t>(data_vec.size());
            write_pod(out, sample_count);
            for(const auto& sample : data_vec) {
                write_pod<T>(out, sample);
            }
        }, signal.data);
    }
}

std::pair<std::unique_ptr<piru::index::VectorSignalStore>, SignalMetadata>
read_signals(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    // --- Common Header ---
    char magic[8];
    in.read(magic, 8);
    if (std::string(magic, 8) != "PIRUSIGS") {
        throw std::runtime_error("Invalid magic number for .signals file.");
    }

    uint32_t format_version;
    read_pod(in, format_version);
    if (format_version > 1000) {
        LOG_WARN("Attempting to read a .signals file with a newer version (" +
                 std::to_string(format_version) + "). Compatibility is not guaranteed.");
    }

    uint32_t header_size;
    read_pod(in, header_size);
    
    std::string backend_type = read_string(in);

    // --- Signal Metadata ---
    SignalMetadata metadata;
    uint64_t node_count;
    read_pod(in, node_count);
    read_pod(in, metadata.quantization_bits);
    read_pod(in, metadata.scale);
    read_pod(in, metadata.offset);
    in.seekg(4, std::ios_base::cur); // Skip reserved

    if (static_cast<uint32_t>(in.tellg()) != header_size) {
        throw std::runtime_error("Header size mismatch in .signals file.");
    }

    auto signals = std::vector<piru::signal::AlignmentQuantizedSignal>();
    signals.reserve(node_count);

    // --- Signal Array ---
    for (uint64_t i = 0; i < node_count; ++i) {
        uint32_t sample_count;
        read_pod(in, sample_count);
        piru::signal::AlignmentQuantizedSignal sig;

        if (backend_type == "float32") {
            sig.kind = piru::signal::AlignmentQuantizationKind::kFloat32;
            std::vector<float> data(sample_count);
            for(uint32_t j=0; j<sample_count; ++j) {
                read_pod(in, data[j]);
            }
            sig.data = std::move(data);
        } else if (backend_type == "int16") {
            sig.kind = piru::signal::AlignmentQuantizationKind::kInt16;
            std::vector<int16_t> data(sample_count);
            for(uint32_t j=0; j<sample_count; ++j) {
                read_pod(in, data[j]);
            }
            sig.data = std::move(data);
        } else if (backend_type == "int8") {
            sig.kind = piru::signal::AlignmentQuantizationKind::kInt8;
            std::vector<int8_t> data(sample_count);
            for(uint32_t j=0; j<sample_count; ++j) {
                read_pod(in, data[j]);
            }
            sig.data = std::move(data);
        } else {
            throw std::runtime_error("Unsupported backend type for .signals file: " + backend_type);
        }
        signals.push_back(std::move(sig));
    }
    
    // Verify that we have reached the end of the file.
    if (in.peek() != EOF) {
        throw std::runtime_error("Trailing data found in .signals file after expected content.");
    }

    auto store = std::make_unique<piru::index::VectorSignalStore>(std::move(signals));
    return {std::move(store), metadata};
}

void write_seeds(const std::string& path, const piru::index::SeedStore& store) {
    const auto* hash_store = dynamic_cast<const piru::index::HashSeedStore*>(&store);
    if (!hash_store) {
        throw std::runtime_error("Unsupported SeedStore backend type for serialization.");
    }

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    // --- Common Header ---
    const char magic[8] = {'P', 'I', 'R', 'U', 'S', 'E', 'E', 'D'};
    out.write(magic, 8);
    write_pod<uint32_t>(out, 1000); // Format version 1.0

    const auto header_size_pos = out.tellp();
    write_pod<uint32_t>(out, 0); // Placeholder for header size
    
    write_string(out, "hash");

    // --- Seed Metadata ---
    const auto& data = hash_store->data();
    uint64_t total_hit_count = 0;
    for (const auto& pair : data) {
        total_hit_count += pair.second.size();
    }

    write_pod<uint64_t>(out, data.size());
    write_pod<uint64_t>(out, total_hit_count);
    write_pod<uint32_t>(out, hash_store->seed_k());
    write_pod<uint32_t>(out, hash_store->seed_stride());
    write_pod<uint32_t>(out, hash_store->seed_qbits());
    write_pod<uint32_t>(out, hash_store->max_hash_frequency());
    write_pod<uint64_t>(out, hash_store->frequency_threshold());
    write_pod<double>(out, hash_store->filter_fraction());

    const uint32_t header_size = static_cast<uint32_t>(out.tellp());
    out.seekp(header_size_pos);
    write_pod<uint32_t>(out, header_size);
    out.seekp(header_size);
    
    // --- Hash Entry Array & Hit List Array ---
    std::vector<uint64_t> sorted_hashes;
    sorted_hashes.reserve(data.size());
    for(const auto& pair : data) {
        sorted_hashes.push_back(pair.first);
    }
    std::sort(sorted_hashes.begin(), sorted_hashes.end());
    
    const auto hash_entry_start = out.tellp();
    const auto hit_list_start = hash_entry_start + std::streampos(sorted_hashes.size() * (sizeof(uint64_t) + sizeof(uint64_t) + sizeof(uint32_t)));
    
    uint64_t current_hit_offset = 0;
    for(const auto& hash : sorted_hashes) {
        const auto& hits = data.at(hash);
        write_pod<uint64_t>(out, hash);
        write_pod<uint64_t>(out, current_hit_offset);
        write_pod<uint32_t>(out, hits.size());
        current_hit_offset += hits.size() * (sizeof(uint32_t) + sizeof(uint32_t));
    }

    out.seekp(hit_list_start);
    for(const auto& hash : sorted_hashes) {
        const auto& hits = data.at(hash);
        for(const auto& hit : hits) {
            write_pod<uint32_t>(out, static_cast<uint32_t>(hit.node_id));
            write_pod<uint32_t>(out, static_cast<uint32_t>(hit.offset));
        }
    }
}

std::unique_ptr<piru::index::HashSeedStore> read_seeds(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    // --- Common Header ---
    char magic[8];
    in.read(magic, 8);
    if (std::string(magic, 8) != "PIRUSEED") {
        throw std::runtime_error("Invalid magic number for .seeds file.");
    }

    uint32_t format_version;
    read_pod(in, format_version);
    if (format_version > 1000) {
        LOG_WARN("Attempting to read a .seeds file with a newer version (" +
                 std::to_string(format_version) + "). Compatibility is not guaranteed.");
    }
    
    uint32_t header_size;
    read_pod(in, header_size);
    std::string backend_type = read_string(in);
    if (backend_type != "hash") {
        throw std::runtime_error("Unsupported backend type for .seeds file: " + backend_type);
    }

    // --- Seed Metadata ---
    uint64_t unique_hash_count, total_hit_count;
    uint32_t seed_k, seed_stride, seed_qbits;
    uint32_t max_hash_frequency;
    uint64_t frequency_threshold;
    double filter_fraction;

    read_pod(in, unique_hash_count);
    read_pod(in, total_hit_count);
    read_pod(in, seed_k);
    read_pod(in, seed_stride);
    read_pod(in, seed_qbits);
    read_pod(in, max_hash_frequency);
    read_pod(in, frequency_threshold);
    read_pod(in, filter_fraction);

    if (static_cast<uint32_t>(in.tellg()) != header_size) {
        throw std::runtime_error("Header size mismatch in .seeds file.");
    }
    
    auto store = std::make_unique<piru::index::HashSeedStore>();
    store->set_seed_k(seed_k);
    store->set_seed_stride(seed_stride);
    store->set_seed_qbits(seed_qbits);
    store->set_max_hash_frequency(max_hash_frequency);
    store->set_frequency_threshold(frequency_threshold);
    store->set_filter_fraction(filter_fraction);
    
    auto& data = store->mutableData();

    struct HashEntry {
        uint64_t hash;
        uint64_t offset;
        uint32_t count;
    };
    std::vector<HashEntry> hash_entries(unique_hash_count);
    for(uint64_t i=0; i < unique_hash_count; ++i) {
        read_pod(in, hash_entries[i].hash);
        read_pod(in, hash_entries[i].offset);
        read_pod(in, hash_entries[i].count);
    }

    const auto hit_list_start = in.tellg();

    for(const auto& entry : hash_entries) {
        in.seekg(hit_list_start + std::streamoff(entry.offset));
        auto& hits = data[entry.hash];
        hits.reserve(entry.count);
        for(uint32_t j=0; j < entry.count; ++j) {
            uint32_t node_id, offset;
            read_pod(in, node_id);
            read_pod(in, offset);
            hits.push_back({node_id, offset});
        }
    }
    
    in.seekg(hit_list_start + std::streamoff(total_hit_count * (sizeof(uint32_t) + sizeof(uint32_t))));
    
    if (in.peek() != EOF) {
        throw std::runtime_error("Trailing data found in .seeds file after expected content.");
    }
    
    return store;
}

} // namespace piru::io::index
