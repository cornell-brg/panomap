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

#include "index/bucket_seed_store.hpp"
#include "index/flat_graph.hpp"
#include "util/logging.hpp"

namespace piru::io::index {

namespace {

constexpr char kMagic[4] = {'P', 'I', 'R', 'X'};
// Version: (major << 16) | minor. Major = breaking format change, minor = compatible.
// NOTE: bump version after first public release. Pre-release format changes
// are breaking but we have no external consumers yet.
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
                const IndexMetadata& metadata,
                const std::vector<float>& node_1d_coords) {
  const auto* adj_store = dynamic_cast<const piru::index::AdjListGraphStore*>(&graph_store);
  if (!adj_store) {
    throw std::runtime_error("Unsupported GraphStore backend for serialization");
  }

  const auto* bucket_store = dynamic_cast<const piru::index::BucketSeedStore*>(&seed_store);
  if (!bucket_store) {
    throw std::runtime_error("Only BucketSeedStore is supported for serialization");
  }

  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open file for writing: " + path);
  }

  const auto& fg = adj_store->flat();

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
  // Legacy: graph_flavor was "vg" vs "dbg", now unused.
  write_string(out, std::string{});
  auto pos_graph = out.tellp();

  /* 3. Graph - nodes (from FlatGraph) */

  write_pod<uint64_t>(out, fg.nodeCount());
  for (std::uint32_t i = 0; i < fg.nodeCount(); ++i) {
    auto name = fg.name(i);
    write_string(out, std::string(name));
    write_pod<uint8_t>(out, fg.isReverse(i) ? 1 : 0);
    if (flags & kFlagHasSequences) {
      write_string(out, fg.seqDecoded(i));
    }
  }

  /* 4. Graph - edges (from FlatGraph CSR) */

  write_pod<uint64_t>(out, fg.edgeCount());
  for (std::uint32_t i = 0; i < fg.nodeCount(); ++i) {
    for (auto it = fg.outBegin(i); it != fg.outEnd(i); ++it) {
      write_pod<uint64_t>(out, i);
      write_pod<uint64_t>(out, *it);
      write_pod<uint64_t>(out, uint64_t{0});  // overlap_bases (not stored in FlatGraph)
    }
  }

  /* 5. Graph - paths (from FlatGraph CSR) */

  write_pod<uint64_t>(out, fg.pathCount());
  for (std::uint32_t i = 0; i < fg.pathCount(); ++i) {
    auto pname = fg.pathName(i);
    write_string(out, std::string(pname));
    write_pod<uint64_t>(out, fg.pathLength(i));
    std::size_t step_count = fg.pathStepCount(i);
    write_pod<uint64_t>(out, step_count);
    const auto* steps = fg.pathStepsBegin(i);
    for (std::size_t j = 0; j < step_count; ++j) {
      write_pod<uint64_t>(out, static_cast<uint64_t>(steps[j]));
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

  /* 7. Seeds (bucket-native format) */

  write_string(out, seed_store.extractor_name());
  const auto& params = seed_store.params();
  write_pod<uint64_t>(out, params.size());
  for (const auto& [key, value] : params) {
    write_string(out, key);
    write_string(out, value);
  }
  write_pod<double>(out, seed_store.filter_fraction());
  write_pod<uint64_t>(out, seed_store.max_hash_frequency());
  write_pod<uint64_t>(out, seed_store.frequency_threshold());

  write_pod<uint32_t>(out, bucket_store->bucket_bits());
  uint64_t num_buckets = bucket_store->num_buckets();
  write_pod<uint64_t>(out, num_buckets);

  for (std::size_t bi = 0; bi < num_buckets; ++bi) {
    const auto& b = bucket_store->bucket(bi);
    uint32_t n_keys = static_cast<uint32_t>(b.keys.size());
    uint32_t n_entries = static_cast<uint32_t>(b.entries.size());
    write_pod<uint32_t>(out, n_keys);
    write_pod<uint32_t>(out, n_entries);
    if (n_keys > 0) {
      out.write(reinterpret_cast<const char*>(b.keys.data()),
                static_cast<std::streamsize>(n_keys * sizeof(uint64_t)));
      out.write(reinterpret_cast<const char*>(b.counts.data()),
                static_cast<std::streamsize>(n_keys * sizeof(uint32_t)));
      out.write(reinterpret_cast<const char*>(b.offsets.data()),
                static_cast<std::streamsize>(n_keys * sizeof(uint32_t)));
    }
    if (n_entries > 0) {
      out.write(reinterpret_cast<const char*>(b.entries.data()),
                static_cast<std::streamsize>(n_entries * sizeof(piru::index::SeedEntry)));
    }
  }

  /* 8. 1D canonical coordinates (optional, float32 bulk write) */

  auto pos_1d = out.tellp();
  uint64_t n_1d = node_1d_coords.size();
  write_pod<uint64_t>(out, n_1d);
  if (n_1d > 0) {
    out.write(reinterpret_cast<const char*>(node_1d_coords.data()),
              static_cast<std::streamsize>(n_1d * sizeof(float)));
  }

  auto pos_end = out.tellp();
  auto total = pos_end - pos_start;
  auto sz_meta = pos_graph - pos_start;
  auto sz_graph = pos_linear - pos_graph;
  auto sz_linear = pos_seeds - pos_linear;
  auto sz_seeds = pos_end - pos_seeds;
  LOG_INFO("Saved index: " + std::to_string(fg.nodeCount()) + " nodes, " +
           std::to_string(seed_store.size()) + " seeds -> " + path);
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

  /* 3-5. Graph - build FlatGraph directly from pirx */

  piru::index::FlatGraph fg;

  // Temporary builder vectors (moved into FlatGraph at the end)
  uint64_t node_count;
  read_pod(in, node_count);

  std::vector<char> seq_data, name_data;
  std::vector<std::uint32_t> seq_offset(node_count), seq_len(node_count);
  std::vector<std::uint32_t> name_offset_nodes(node_count);
  std::vector<std::uint16_t> name_len_nodes(node_count);
  std::vector<std::uint8_t> is_reverse(node_count);

  for (uint64_t i = 0; i < node_count; ++i) {
    std::string orig_id = read_string(in);
    uint8_t rev;
    read_pod(in, rev);
    is_reverse[i] = rev;

    name_offset_nodes[i] = static_cast<std::uint32_t>(name_data.size());
    name_len_nodes[i] = static_cast<std::uint16_t>(orig_id.size());
    name_data.insert(name_data.end(), orig_id.begin(), orig_id.end());

    if (flags & kFlagHasSequences) {
      std::string seq = read_string(in);
      seq_offset[i] = static_cast<std::uint32_t>(seq_data.size());
      seq_len[i] = static_cast<std::uint32_t>(seq.size());
      seq_data.insert(seq_data.end(), seq.begin(), seq.end());
    }
  }

  /* 4. Edges */
  uint64_t edge_count;
  read_pod(in, edge_count);

  // Build CSR from edge list
  std::vector<std::vector<std::uint32_t>> adj(node_count);
  for (uint64_t i = 0; i < edge_count; ++i) {
    uint64_t from, to, overlap;
    read_pod(in, from);
    read_pod(in, to);
    read_pod(in, overlap);
    if (from < node_count) {
      adj[from].push_back(static_cast<std::uint32_t>(to));
    }
  }

  std::vector<std::uint32_t> edge_target;
  std::vector<std::uint32_t> out_edge_offset(node_count + 1);
  for (uint64_t i = 0; i < node_count; ++i) {
    out_edge_offset[i] = static_cast<std::uint32_t>(edge_target.size());
    edge_target.insert(edge_target.end(), adj[i].begin(), adj[i].end());
  }
  out_edge_offset[node_count] = static_cast<std::uint32_t>(edge_target.size());
  adj.clear();

  /* 5. Paths */
  uint64_t path_count;
  read_pod(in, path_count);

  std::vector<std::uint32_t> step_data;
  std::vector<std::uint32_t> path_step_offset(path_count + 1);
  std::vector<std::uint32_t> path_name_offset(path_count);
  std::vector<std::uint16_t> path_name_len(path_count);
  std::vector<std::uint64_t> path_length(path_count);

  for (uint64_t i = 0; i < path_count; ++i) {
    std::string pname = read_string(in);
    path_name_offset[i] = static_cast<std::uint32_t>(name_data.size());
    path_name_len[i] = static_cast<std::uint16_t>(pname.size());
    name_data.insert(name_data.end(), pname.begin(), pname.end());

    read_pod(in, path_length[i]);

    uint64_t step_count;
    read_pod(in, step_count);
    path_step_offset[i] = static_cast<std::uint32_t>(step_data.size());
    for (uint64_t j = 0; j < step_count; ++j) {
      uint64_t nid;
      read_pod(in, nid);
      step_data.push_back(static_cast<std::uint32_t>(nid));
    }
  }
  path_step_offset[path_count] = static_cast<std::uint32_t>(step_data.size());

  // Assemble FlatGraph
  fg = piru::index::FlatGraph::fromRawArrays(
      static_cast<std::uint32_t>(node_count), static_cast<std::uint32_t>(path_count),
      std::move(seq_data), std::move(seq_offset), std::move(seq_len),
      std::move(name_data), std::move(name_offset_nodes), std::move(name_len_nodes),
      std::move(is_reverse),
      std::move(edge_target), std::move(out_edge_offset),
      std::move(step_data), std::move(path_step_offset),
      std::move(path_name_offset), std::move(path_name_len), std::move(path_length));

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

  /* 7. Seeds (bucket-native format) */

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

  uint32_t bucket_bits;
  read_pod(in, bucket_bits);

  uint64_t num_buckets;
  read_pod(in, num_buckets);

  std::vector<piru::index::Bucket> buckets(num_buckets);
  std::size_t total_hits = 0;

  for (uint64_t bi = 0; bi < num_buckets; ++bi) {
    uint32_t n_keys, n_entries;
    read_pod(in, n_keys);
    read_pod(in, n_entries);

    auto& b = buckets[bi];
    if (n_keys > 0) {
      b.keys.resize(n_keys);
      b.counts.resize(n_keys);
      b.offsets.resize(n_keys);
      in.read(reinterpret_cast<char*>(b.keys.data()),
              static_cast<std::streamsize>(n_keys * sizeof(uint64_t)));
      in.read(reinterpret_cast<char*>(b.counts.data()),
              static_cast<std::streamsize>(n_keys * sizeof(uint32_t)));
      in.read(reinterpret_cast<char*>(b.offsets.data()),
              static_cast<std::streamsize>(n_keys * sizeof(uint32_t)));
    }
    if (n_entries > 0) {
      b.entries.resize(n_entries);
      in.read(reinterpret_cast<char*>(b.entries.data()),
              static_cast<std::streamsize>(n_entries * sizeof(piru::index::SeedEntry)));
    }
    total_hits += n_entries;
  }

  auto seeds = std::make_unique<piru::index::BucketSeedStore>(
      std::move(buckets), bucket_bits,
      std::move(seed_extractor_name), std::move(seed_params),
      max_freq, freq_threshold, filter_fraction);

  /* 8. 1D canonical coordinates (optional, backwards-compatible) */

  std::vector<float> node_1d_coords;
  if (in.peek() != std::ifstream::traits_type::eof()) {
    uint64_t n_1d;
    read_pod(in, n_1d);
    if (n_1d > 0) {
      node_1d_coords.resize(n_1d);
      in.read(reinterpret_cast<char*>(node_1d_coords.data()),
              static_cast<std::streamsize>(n_1d * sizeof(float)));
      LOG_INFO("Loaded 1D coords: " + std::to_string(n_1d) + " nodes");
    }
  }

  LOG_INFO("Loaded index: " + std::to_string(node_count) + " nodes, " +
           std::to_string(seeds->size()) + " seeds (" + std::to_string(total_hits) +
           " hits) ← " + path);

  return {std::move(metadata),
          std::make_unique<piru::index::FlatGraphStore>(std::move(fg)),
          std::move(seeds), std::move(linearization_coords), std::move(node_1d_coords)};
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
