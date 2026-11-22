/**
 * inspect.cpp
 *
 * CLI handler for `panomap inspect`. Loads a .pirx index and prints
 * metadata, graph stats, and seed store stats.
 *
 * SPDX-License-Identifier: MIT
 */

#include "signal/commands/inspect.hpp"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "core/io/index/serialization.hpp"
#include "core/util/logging.hpp"

namespace {

// True if bytes [off, off + len) are the same in both files.
bool same_bytes(const std::string& a, const std::string& b, std::uint64_t off, std::uint64_t len) {
  std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
  fa.seekg(static_cast<std::streamoff>(off));
  fb.seekg(static_cast<std::streamoff>(off));
  std::vector<char> ba(1 << 20), bb(1 << 20);
  while (len > 0) {
    const auto n = static_cast<std::streamsize>(std::min<std::uint64_t>(len, ba.size()));
    fa.read(ba.data(), n);
    fb.read(bb.data(), n);
    if (fa.gcount() != n || fb.gcount() != n) return false;
    if (std::memcmp(ba.data(), bb.data(), static_cast<std::size_t>(n)) != 0) return false;
    len -= static_cast<std::uint64_t>(n);
  }
  return true;
}

// Prints identical or differs for each section of two indexes. Two 32-thread builds of
// the same graph differ only in coords_1d: the 1D layout runs lock-free parallel SGD
// (Hogwild, as odgi's PG-SGD), so its coordinates change slightly from run to run.
int compare_indexes(const std::string& a, const panomap::io::index::LoadedIndex& la,
                    const std::string& b) {
  auto lb = panomap::io::index::load_index(b);
  const auto& ma = la.metadata;
  const auto& mb = lb.metadata;
  const auto& sa = la.section_sizes;
  const auto& sb = lb.section_sizes;
  int differ = 0;
  auto report = [&](const char* name, bool same, const char* note = "") {
    std::cout << "  " << std::left << std::setw(13) << name << (same ? "identical" : "differs")
              << note << "\n";
    differ += same ? 0 : 1;
  };
  std::cout << "compare: " << a << "\n     vs: " << b << "\n";
  report("header_meta",
         ma.mode == mb.mode && ma.version == mb.version && ma.model_name == mb.model_name &&
             ma.pore_k == mb.pore_k && ma.tokenizer == mb.tokenizer &&
             sa.header_meta == sb.header_meta,
         ma.build_timestamp != mb.build_timestamp ? " (build time not compared)" : "");
  struct Sec { const char* name; std::uint64_t size_a, size_b; };
  const Sec secs[] = {{"graph", sa.graph, sb.graph},
                      {"path_steps", sa.path_steps, sb.path_steps},
                      {"seeds", sa.seeds, sb.seeds},
                      {"coords_1d", sa.coords_1d, sb.coords_1d},
                      {"components", sa.components, sb.components}};
  std::uint64_t off = sa.header_meta;
  for (const auto& sec : secs) {
    report(sec.name, sec.size_a == sec.size_b && sa.header_meta == sb.header_meta &&
                         same_bytes(a, b, off, sec.size_a));
    off += sec.size_a;
  }
  return differ == 0 ? 0 : 2;
}

}  // namespace

int handle_inspect(const std::vector<std::string>& args) {
  panomap::cli::Parsed parsed;
  panomap::cli::ParseConfig config;
  config.usage = "Usage: panomap inspect [--compare other.pirx] <index.pirx>";
  config.positional_help = {"<index.pirx>    Index file to inspect"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'\0', "dump-path-coords", true, "Dump per-path canonical coordinate mapping to TSV"},
      {'\0', "compare", true,
       "Compare with another index, section by section, and print identical or differs"},
  };
  config.on_error = [](const std::string&) { std::cerr << "inspect: invalid option\n"; };

  if (!panomap::cli::parse_args(args, config, parsed)) {
    panomap::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    panomap::cli::print_help(config, std::cout);
    return 0;
  }
  if (parsed.positionals.empty()) {
    LOG_ERROR("inspect: missing required <index.pirx>");
    panomap::cli::print_help(config, std::cerr);
    return 1;
  }

  const std::string index_path = parsed.positionals[0];

  if (!panomap::io::index::is_pirx_index(index_path)) {
    LOG_ERROR("inspect: not a valid .pirx file: " + index_path);
    return 1;
  }

  auto loaded = panomap::io::index::load_index(index_path);

  if (loaded.metadata.mode != panomap::io::index::IndexMode::kSignal) {
    LOG_ERROR(std::string("inspect: index was built in mode '") +
              panomap::io::index::mode_name(loaded.metadata.mode) +
              "', but panomap only loads 'signal' indexes. Use panomap-base instead.");
    return 1;
  }

  if (parsed.values.count("compare")) {
    return compare_indexes(index_path, loaded, parsed.values.at("compare"));
  }

  const auto& meta = loaded.metadata;
  const auto& fg = loaded.graph->flat();
  const auto& seeds = *loaded.seeds;

  const auto& params = seeds.params();
  auto get = [&](const std::string& key) -> std::string {
    auto it = params.find(key);
    return it != params.end() ? it->second : "";
  };

  std::cout << "index: " << index_path << "\n";
  std::cout << "panomap_version: " << (meta.version.empty() ? "unknown" : meta.version) << "\n";
  if (meta.build_timestamp > 0) {
    auto t = static_cast<std::time_t>(meta.build_timestamp);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::cout << "built: " << buf << "\n";
  }
  std::cout << "model: " << meta.model_name << "\n";
  std::cout << "pore_k: " << meta.pore_k << "\n";
  std::cout << "tokenizer: " << meta.tokenizer << "\n";
  std::cout << "seed_extractor: " << seeds.extractor_name() << "\n";
  std::cout << "seed_k: " << get("k") << "\n";
  std::cout << "seed_qbits: " << get("qbits") << "\n";
  std::cout << "token_collapse: " << get("collapse") << "\n";
  std::cout << "diff_filter: " << get("diff") << "\n";
  std::cout << "preset: " << get("preset") << "\n";
  std::cout << "\n";
  std::cout << "graph:\n";
  std::cout << "  nodes: " << fg.nodeCount() << "\n";
  std::cout << "  edges: " << fg.edgeCount() << "\n";
  std::cout << "  paths:\n";
  for (std::uint32_t i = 0; i < fg.pathCount(); ++i) {
    std::cout << "    - " << fg.pathName(i) << " (" << fg.pathLength(i) << " bp)\n";
  }
  std::cout << "\n";
  std::cout << "seeds:\n";
  std::cout << "  unique_hashes: " << seeds.size() << "\n";
  std::cout << "  max_frequency: " << seeds.max_hash_frequency() << "\n";
  std::cout << "  freq_threshold: " << seeds.frequency_threshold() << "\n";
  std::cout << "\n";
  std::cout << "path_steps: " << fg.pathStepOffsets().back() << "\n";
  std::cout << "1d_coordinates: " << (loaded.node_1d_coords.empty() ? "no" : "yes") << "\n";
  std::cout << "component_ids: " << (loaded.component_ids.empty() ? "no" : "yes") << "\n";

  /* On-disk section breakdown (bytes + % of total). */
  const auto& sz = loaded.section_sizes;
  auto pct = [&](uint64_t v) {
    return sz.total > 0 ? (100.0 * static_cast<double>(v) / static_cast<double>(sz.total)) : 0.0;
  };
  std::cout << "\n";
  std::cout << "sections (bytes, % of total):\n";
  std::cout << "  total:         " << sz.total << "\n";
  std::cout << "  header_meta:   " << sz.header_meta << "  (" << std::fixed << std::setprecision(1)
            << pct(sz.header_meta) << "%)\n";
  std::cout << "  graph:         " << sz.graph << "  (" << pct(sz.graph) << "%)\n";
  std::cout << "  path_steps:    " << sz.path_steps << "  (" << pct(sz.path_steps) << "%)\n";
  std::cout << "  seeds:         " << sz.seeds << "  (" << pct(sz.seeds) << "%)\n";
  std::cout << "  coords_1d:     " << sz.coords_1d << "  (" << pct(sz.coords_1d) << "%)\n";
  std::cout << "  components:    " << sz.components << "  (" << pct(sz.components) << "%)\n";

  /* Dump per-node canonical coordinates.
   * One row per node: node_id, canon_start (fwd), canon_end (rev), component_id.
   * Eval scripts use this + GFA paths to map base-space positions to canonical space. */
  if (parsed.values.count("dump-path-coords")) {
    if (loaded.node_1d_coords.empty() || loaded.component_ids.empty()) {
      LOG_ERROR("dump-path-coords requires 1D coordinates and component IDs in the index");
      return 1;
    }
    const auto& coords_1d = loaded.node_1d_coords;
    const auto& comp_ids = loaded.component_ids;

    std::ofstream ofs(parsed.values.at("dump-path-coords"));
    if (!ofs) {
      LOG_ERROR("Failed to open output: " + parsed.values.at("dump-path-coords"));
      return 1;
    }
    ofs << "node_id\tcanon_start\tcanon_end\tcomponent_id\n";
    for (std::uint32_t i = 0; i < fg.nodeCount(); i += 2) {
      ofs << i << '\t' << coords_1d[i] << '\t' << coords_1d[i + 1] << '\t' << comp_ids[i] << '\n';
    }
    LOG_INFO("Dumped node coords: " + std::to_string(fg.nodeCount() / 2) + " nodes to " +
             parsed.values.at("dump-path-coords"));
  }

  return 0;
}
