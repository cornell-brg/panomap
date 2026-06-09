/**
 * serialization.cpp
 *
 * Binary read/write for the .pirx index format.
 * See serialization.hpp for format layout.
 *
 * SPDX-License-Identifier: MIT
 */

#include "core/io/index/serialization.hpp"

#include <fstream>
#include <stdexcept>

#include "core/index/bucket_seed_store.hpp"
#include "core/index/flat_graph.hpp"
#include "core/util/logging.hpp"

namespace panomap::io::index {

namespace {

constexpr char kMagic[4] = {'P', 'I', 'R', 'X'};
// Version: (major << 16) | minor. Major = breaking format change, minor = compatible.
// v2: added mode byte after flags (signal vs. base).
// v2.1 (dev-112 Phase 1): drop isReverse byte from per-node record. Derived
// from canonical pair convention (fwd=2i, rev=2i+1) instead.
// v2.2 (dev-112 Phase 2): drop edges section. Chainer never iterates edges
// at map time; sort_1d / gfa_exporter only see edges in the build-time
// in-memory graph. Loaded indexes get an empty CSR (out_edge_offset[N+1]=0).
// v2.3 (dev-112 Phase 3): drop path steps array. Map-time consumers
// (gaf_writer) only walk steps as a fallback when pathLength is 0 -- which
// can't happen for indexes built through the standard pipeline. Loaded
// indexes get empty step_data + zeroed path_step_offset.
// Minor-bump convention: any incremental serialization change inside the v2.x
// series increments kVersionMinor. The load-time check requires strict full
// version equality, so older v2.x indexes are refused (force re-index).
// Reserve major bumps for bigger semantic breaks.
constexpr uint32_t kVersionMajor = 2;
constexpr uint32_t kVersionMinor = 3;
constexpr uint32_t kVersion = (kVersionMajor << 16) | kVersionMinor;

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

const char* mode_name(IndexMode mode) {
  switch (mode) {
    case IndexMode::kSignal: return "signal";
    case IndexMode::kBase:   return "base";
  }
  return "unknown";
}

void save_index(const std::string& path, const panomap::index::GraphStore& graph_store,
                const panomap::index::SeedStore& seed_store,
                const std::vector<std::vector<panomap::index::LinearCoordinate>>& linearization_coords,
                const IndexMetadata& metadata, const std::vector<float>& node_1d_coords,
                const std::vector<std::uint32_t>& component_ids) {
  const auto* adj_store = dynamic_cast<const panomap::index::AdjListGraphStore*>(&graph_store);
  if (!adj_store) {
    throw std::runtime_error("Unsupported GraphStore backend for serialization");
  }

  const auto* bucket_store = dynamic_cast<const panomap::index::BucketSeedStore*>(&seed_store);
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
  write_pod<uint8_t>(out, static_cast<uint8_t>(metadata.mode));

  /* 2. Metadata */

  write_string(out, metadata.version);
  write_pod<uint64_t>(out, metadata.build_timestamp);
  write_string(out, metadata.model_name);
  write_pod<uint32_t>(out, metadata.pore_k);
  write_string(out, metadata.tokenizer);
  auto pos_graph = out.tellp();

  /* 3. Graph - nodes (from FlatGraph) */

  // dev-112: isReverse not serialized; derived from (node_id & 1) at load time.
  // Canonical pair convention enforced by simple_expand (fwd=2i, rev=2i+1).
  write_pod<uint64_t>(out, fg.nodeCount());
  for (std::uint32_t i = 0; i < fg.nodeCount(); ++i) {
    auto name = fg.name(i);
    write_string(out, std::string(name));
    write_pod<uint32_t>(out, static_cast<uint32_t>(fg.seqLen(i)));
  }

  /* 4. Graph - edges: not serialized as of v2.2 (dev-112 Phase 2).
   * The chainer never iterates edges at map time; sort_1d / gfa_exporter
   * only need them at index build (in-memory). Dropping them saves
   * ~24 bytes per edge (yeast N=8: ~42 MB, zymo N=8: ~183 MB). */

  /* 5. Graph - paths: name + length only (v2.3, dev-112 Phase 3).
   * Path steps array dropped from disk -- only gfa_exporter (build-time
   * only) and the gaf_writer length fallback (never fires when pathLength
   * is populated by the indexer pipeline, which is the only path that
   * writes pirx files) consume them. */

  write_pod<uint64_t>(out, fg.pathCount());
  for (std::uint32_t i = 0; i < fg.pathCount(); ++i) {
    auto pname = fg.pathName(i);
    write_string(out, std::string(pname));
    write_pod<uint64_t>(out, fg.pathLength(i));
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
                static_cast<std::streamsize>(n_entries * sizeof(panomap::index::SeedEntry)));
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

  /* 9. Connected component IDs (optional, uint32 bulk write) */

  uint64_t n_comp = component_ids.size();
  write_pod<uint64_t>(out, n_comp);
  if (n_comp > 0) {
    out.write(reinterpret_cast<const char*>(component_ids.data()),
              static_cast<std::streamsize>(n_comp * sizeof(uint32_t)));
  }

  auto pos_end = out.tellp();
  auto total = pos_end - pos_start;
  auto sz_meta = pos_graph - pos_start;
  auto sz_graph = pos_linear - pos_graph;
  auto sz_linear = pos_seeds - pos_linear;
  auto sz_seeds = pos_end - pos_seeds;
  LOG_INFO("Saved index: " + std::to_string(fg.nodeCount()) + " nodes, " +
           std::to_string(seed_store.size()) + " seeds -> " + path);
  LOG_INFO("Index breakdown: total=" + std::to_string(total / (1024 * 1024)) + "MB" +
           "  header+meta=" + std::to_string(sz_meta * 100 / total) + "%" +
           "  graph=" + std::to_string(sz_graph * 100 / total) + "%" + " (" +
           std::to_string(sz_graph / (1024 * 1024)) + "MB)" +
           "  linearization=" + std::to_string(sz_linear * 100 / total) + "%" + " (" +
           std::to_string(sz_linear / (1024 * 1024)) + "MB)" +
           "  seeds=" + std::to_string(sz_seeds * 100 / total) + "%" + " (" +
           std::to_string(sz_seeds / (1024 * 1024)) + "MB)");
}

LoadedIndex load_index(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open file for reading: " + path);
  }

  // Track byte positions at each section boundary (matches save_index layout).
  std::int64_t pos_start = 0, pos_graph = 0, pos_linear = 0, pos_seeds = 0,
               pos_1d = 0, pos_components = 0, pos_end = 0;

  /* 1. Header */

  char magic[4];
  in.read(magic, 4);
  if (std::string(magic, 4) != std::string(kMagic, 4)) {
    throw std::runtime_error("Invalid magic number for index file");
  }

  uint32_t version;
  read_pod(in, version);
  if (version != kVersion) {
    uint32_t major = version >> 16;
    uint32_t minor = version & 0xFFFF;
    throw std::runtime_error("Incompatible index version " + std::to_string(major) + "." +
                             std::to_string(minor) + " (expected " +
                             std::to_string(kVersionMajor) + "." +
                             std::to_string(kVersionMinor) + ", re-index required)");
  }

  uint32_t flags;
  read_pod(in, flags);

  uint8_t mode_byte;
  read_pod(in, mode_byte);
  if (mode_byte > static_cast<uint8_t>(IndexMode::kBase)) {
    throw std::runtime_error("Unknown index mode byte: " + std::to_string(mode_byte));
  }

  /* 2. Metadata */

  IndexMetadata metadata;
  metadata.mode = static_cast<IndexMode>(mode_byte);
  metadata.version = read_string(in);
  read_pod(in, metadata.build_timestamp);
  metadata.model_name = read_string(in);
  read_pod(in, metadata.pore_k);
  metadata.tokenizer = read_string(in);

  pos_graph = static_cast<std::int64_t>(in.tellg());

  /* 3-5. Graph - build FlatGraph directly from pirx */

  panomap::index::FlatGraph fg;

  // Temporary builder vectors (moved into FlatGraph at the end)
  uint64_t node_count;
  read_pod(in, node_count);

  std::vector<char> seq_data, name_data;
  std::vector<std::uint32_t> seq_offset(node_count), seq_len(node_count);
  std::vector<std::uint32_t> name_offset_nodes(node_count);
  std::vector<std::uint16_t> name_len_nodes(node_count);

  // Phase 1 (dev-112): isReverse is derived from (node_id & 1) instead of a
  // stored byte. Canonical pair convention from simple_expand: fwd=2i, rev=2i+1.
  std::vector<std::uint32_t> saved_seq_len(node_count);
  for (uint64_t i = 0; i < node_count; ++i) {
    std::string orig_id = read_string(in);
    uint32_t slen;
    read_pod(in, slen);
    saved_seq_len[i] = slen;

    name_offset_nodes[i] = static_cast<std::uint32_t>(name_data.size());
    name_len_nodes[i] = static_cast<std::uint16_t>(orig_id.size());
    name_data.insert(name_data.end(), orig_id.begin(), orig_id.end());
  }

  /* 4. Edges: not serialized in v2.2+ (dev-112 Phase 2). The chainer doesn't
   * traverse edges at map time, so we create an empty CSR. outBegin/outEnd
   * return zero-length ranges; outDegree returns 0. Build-time consumers
   * (sort_1d / gfa_exporter) only see edges when constructing the in-memory
   * graph -- they never come from a loaded index. */
  std::vector<std::uint32_t> edge_target;
  std::vector<std::uint32_t> out_edge_offset(node_count + 1, 0);

  /* 5. Paths: name + length only in v2.3+ (dev-112 Phase 3). Steps array
   * dropped from disk; we synthesize empty offsets so the accessors keep
   * working (pathStepsBegin/End return zero-length ranges). The gaf_writer
   * length fallback uses these only when pathLength is 0 -- which can't
   * happen for indexes built through the standard pipeline. */
  uint64_t path_count;
  read_pod(in, path_count);

  std::vector<std::uint32_t> step_data;  // empty
  std::vector<std::uint32_t> path_step_offset(path_count + 1, 0);
  std::vector<std::uint32_t> path_name_offset(path_count);
  std::vector<std::uint16_t> path_name_len(path_count);
  std::vector<std::uint64_t> path_length(path_count);

  for (uint64_t i = 0; i < path_count; ++i) {
    std::string pname = read_string(in);
    path_name_offset[i] = static_cast<std::uint32_t>(name_data.size());
    path_name_len[i] = static_cast<std::uint16_t>(pname.size());
    name_data.insert(name_data.end(), pname.begin(), pname.end());

    read_pod(in, path_length[i]);
  }

  // Assemble FlatGraph (seq_len passed as zeros; real lengths set below)
  fg = panomap::index::FlatGraph::fromRawArrays(
      static_cast<std::uint32_t>(node_count), static_cast<std::uint32_t>(path_count),
      std::move(seq_data), std::move(seq_offset), std::move(seq_len), std::move(name_data),
      std::move(name_offset_nodes), std::move(name_len_nodes),
      std::move(edge_target), std::move(out_edge_offset), std::move(step_data),
      std::move(path_step_offset), std::move(path_name_offset), std::move(path_name_len),
      std::move(path_length));
  fg.setSeqLens(std::move(saved_seq_len));

  pos_linear = static_cast<std::int64_t>(in.tellg());

  /* 6. Linearization */

  uint64_t lin_node_count;
  read_pod(in, lin_node_count);

  std::vector<std::vector<panomap::index::LinearCoordinate>> linearization_coords(lin_node_count);
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

  pos_seeds = static_cast<std::int64_t>(in.tellg());

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

  std::vector<panomap::index::Bucket> buckets(num_buckets);
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
              static_cast<std::streamsize>(n_entries * sizeof(panomap::index::SeedEntry)));
    }
    total_hits += n_entries;
  }

  auto seeds = std::make_unique<panomap::index::BucketSeedStore>(
      std::move(buckets), bucket_bits, std::move(seed_extractor_name), std::move(seed_params),
      max_freq, freq_threshold, filter_fraction);

  pos_1d = static_cast<std::int64_t>(in.tellg());

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

  pos_components = static_cast<std::int64_t>(in.tellg());

  /* 9. Connected component IDs (optional, backwards-compatible) */

  std::vector<std::uint32_t> component_ids;
  if (in.peek() != std::ifstream::traits_type::eof()) {
    uint64_t n_comp;
    read_pod(in, n_comp);
    if (n_comp > 0) {
      component_ids.resize(n_comp);
      in.read(reinterpret_cast<char*>(component_ids.data()),
              static_cast<std::streamsize>(n_comp * sizeof(uint32_t)));
      LOG_INFO("Loaded component IDs: " + std::to_string(n_comp) + " nodes");
    }
  }

  LOG_INFO("Loaded index: " + std::to_string(node_count) + " nodes, " +
           std::to_string(seeds->size()) + " seeds (" + std::to_string(total_hits) + " hits) ← " +
           path);

  // tellg() returns -1 at EOF; resolve to actual file size for the final section.
  in.clear();
  in.seekg(0, std::ios::end);
  pos_end = static_cast<std::int64_t>(in.tellg());

  SectionSizes sz;
  sz.header_meta = static_cast<uint64_t>(pos_graph - pos_start);
  sz.graph = static_cast<uint64_t>(pos_linear - pos_graph);
  sz.linearization = static_cast<uint64_t>(pos_seeds - pos_linear);
  sz.seeds = static_cast<uint64_t>(pos_1d - pos_seeds);
  sz.coords_1d = static_cast<uint64_t>(pos_components - pos_1d);
  sz.components = static_cast<uint64_t>(pos_end - pos_components);
  sz.total = static_cast<uint64_t>(pos_end - pos_start);

  return {std::move(metadata), std::make_unique<panomap::index::FlatGraphStore>(std::move(fg)),
          std::move(seeds), std::move(linearization_coords), std::move(node_1d_coords),
          std::move(component_ids), sz};
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

}  // namespace panomap::io::index
