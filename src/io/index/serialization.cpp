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
// Version: (major << 16) | minor. Major = breaking format change, minor = compatible.
constexpr uint32_t kVersionMajor = 1;
constexpr uint32_t kVersionMinor = 0;
constexpr uint32_t kVersion = (kVersionMajor << 16) | kVersionMinor;

constexpr uint32_t kFlagHasSequences = 1 << 0;
constexpr uint32_t kFlagHasFuzzySignals = 1 << 1;  // reserved
constexpr uint32_t kFlagHasAlignSignals = 1 << 2;  // reserved

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

  auto pos_start = out.tellp();
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
  auto pos_graph = out.tellp();

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

  auto pos_linear = out.tellp();

  /* 6. Linearization */

  write_pod<uint64_t>(out, linearization_coords.size());
  for (const auto& coords : linearization_coords) {
    write_pod<uint64_t>(out, coords.size());
    for (const auto& coord : coords) {
      write_pod<uint64_t>(out, coord.path_id);
      write_pod<int64_t>(out, coord.ref_coord);
    }
  }

  auto pos_seeds = out.tellp();

  /* 7. Seeds (CSR format) */

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

  // Build CSR from HashSeedStore: sort hashes, flatten hits
  const auto& seed_data = hash_store->data();
  std::vector<std::uint64_t> sorted_hashes;
  sorted_hashes.reserve(seed_data.size());
  for (const auto& [hash, hits] : seed_data) {
    sorted_hashes.push_back(hash);
  }
  std::sort(sorted_hashes.begin(), sorted_hashes.end());

  std::vector<std::uint32_t> csr_offsets;
  csr_offsets.reserve(sorted_hashes.size() + 1);
  std::vector<piru::index::SeedEntry> csr_entries;

  std::uint32_t running_offset = 0;
  for (auto h : sorted_hashes) {
    csr_offsets.push_back(running_offset);
    const auto& hits = seed_data.at(h);
    csr_entries.insert(csr_entries.end(), hits.begin(), hits.end());
    running_offset += static_cast<std::uint32_t>(hits.size());
  }
  csr_offsets.push_back(running_offset);  // sentinel

  uint64_t n_hashes = sorted_hashes.size();
  uint64_t n_entries = csr_entries.size();
  write_pod<uint64_t>(out, n_hashes);
  write_pod<uint64_t>(out, n_entries);
  out.write(reinterpret_cast<const char*>(sorted_hashes.data()),
            static_cast<std::streamsize>(n_hashes * sizeof(uint64_t)));
  out.write(reinterpret_cast<const char*>(csr_offsets.data()),
            static_cast<std::streamsize>((n_hashes + 1) * sizeof(uint32_t)));
  out.write(reinterpret_cast<const char*>(csr_entries.data()),
            static_cast<std::streamsize>(n_entries * sizeof(piru::index::SeedEntry)));

  auto pos_end = out.tellp();
  auto total = pos_end - pos_start;
  auto sz_meta = pos_graph - pos_start;
  auto sz_graph = pos_linear - pos_graph;
  auto sz_linear = pos_seeds - pos_linear;
  auto sz_seeds = pos_end - pos_seeds;
  LOG_INFO("Saved index: " + std::to_string(graph.nodeCount()) + " nodes, " +
           std::to_string(n_hashes) + " seeds (" + std::to_string(n_entries) + " hits) -> " + path);
  LOG_INFO("Index breakdown: total=" + std::to_string(total / (1024*1024)) + "MB" +
           "  header+meta=" + std::to_string(sz_meta * 100 / total) + "%" +
           "  graph=" + std::to_string(sz_graph * 100 / total) + "%" +
           " (" + std::to_string(sz_graph / (1024*1024)) + "MB)" +
           "  linearization=" + std::to_string(sz_linear * 100 / total) + "%" +
           " (" + std::to_string(sz_linear / (1024*1024)) + "MB)" +
           "  seeds=" + std::to_string(sz_seeds * 100 / total) + "%" +
           " (" + std::to_string(sz_seeds / (1024*1024)) + "MB)");
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
  uint32_t major = version >> 16;
  if (major != kVersionMajor) {
    throw std::runtime_error("Incompatible index version " +
        std::to_string(major) + "." + std::to_string(version & 0xFFFF) +
        " (expected " + std::to_string(kVersionMajor) + ".x)");
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

  /* 7. Seeds (CSR format) */

  std::string seed_extractor_name = read_string(in);

  uint64_t param_count;
  read_pod(in, param_count);
  std::map<std::string, std::string> seed_params;
  for (uint64_t i = 0; i < param_count; ++i) {
    std::string key = read_string(in);
    std::string value = read_string(in);
    seed_params[key] = value;
  }

  double filter_fraction;
  read_pod(in, filter_fraction);

  uint64_t max_freq;
  read_pod(in, max_freq);

  uint64_t freq_threshold;
  read_pod(in, freq_threshold);

  // Bulk-read CSR arrays
  uint64_t n_hashes, n_entries;
  read_pod(in, n_hashes);
  read_pod(in, n_entries);

  std::vector<uint64_t> hashes(n_hashes);
  in.read(reinterpret_cast<char*>(hashes.data()),
          static_cast<std::streamsize>(n_hashes * sizeof(uint64_t)));

  std::vector<uint32_t> offsets(n_hashes + 1);
  in.read(reinterpret_cast<char*>(offsets.data()),
          static_cast<std::streamsize>((n_hashes + 1) * sizeof(uint32_t)));

  std::vector<piru::index::SeedEntry> entries(n_entries);
  in.read(reinterpret_cast<char*>(entries.data()),
          static_cast<std::streamsize>(n_entries * sizeof(piru::index::SeedEntry)));

  auto seeds = std::make_unique<piru::index::FlatSeedStore>(
      std::move(hashes), std::move(offsets), std::move(entries),
      std::move(seed_extractor_name), std::move(seed_params),
      max_freq, freq_threshold, filter_fraction);

  LOG_INFO("Loaded index: " + std::to_string(node_count) + " nodes, " +
           std::to_string(n_hashes) + " seeds (" + std::to_string(n_entries) + " hits) ← " + path);

  return {std::move(metadata), std::make_unique<piru::index::AdjListGraphStore>(std::move(*graph)),
          std::move(seeds), std::move(linearization_coords)};
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
