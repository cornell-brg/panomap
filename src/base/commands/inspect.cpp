/**
 * inspect.cpp (base mode)
 *
 * CLI handler for `piru-base inspect`. Loads a .pirx index, refuses
 * non-base indexes, and prints metadata + graph + seed stats.
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/commands/inspect.hpp"

#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "core/io/index/serialization.hpp"
#include "core/util/logging.hpp"

int handle_base_inspect(const std::vector<std::string>& args) {
  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru-base inspect <index.pirx>";
  config.positional_help = {"<index.pirx>    Index file to inspect"};
  config.options = {
      {'h', "help", false, "Show help"},
      {'\0', "dump-path-coords", true, "Dump per-node canonical coordinate mapping to TSV"},
  };
  config.on_error = [](const std::string&) { std::cerr << "inspect: invalid option\n"; };

  if (!piru::cli::parse_args(args, config, parsed)) {
    piru::cli::print_help(config, std::cerr);
    return 1;
  }
  if (parsed.values.count("help")) {
    piru::cli::print_help(config, std::cout);
    return 0;
  }
  if (parsed.positionals.empty()) {
    LOG_ERROR("inspect: missing required <index.pirx>");
    piru::cli::print_help(config, std::cerr);
    return 1;
  }

  const std::string index_path = parsed.positionals[0];

  if (!piru::io::index::is_pirx_index(index_path)) {
    LOG_ERROR("inspect: not a valid .pirx file: " + index_path);
    return 1;
  }

  auto loaded = piru::io::index::load_index(index_path);

  if (loaded.metadata.mode != piru::io::index::IndexMode::kBase) {
    LOG_ERROR(std::string("inspect: index was built in mode '") +
              piru::io::index::mode_name(loaded.metadata.mode) +
              "', but piru-base only loads 'base' indexes. Use piru-signal instead.");
    return 1;
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
  std::cout << "mode: " << piru::io::index::mode_name(meta.mode) << "\n";
  std::cout << "piru_version: " << (meta.version.empty() ? "unknown" : meta.version) << "\n";
  if (meta.build_timestamp > 0) {
    auto t = static_cast<std::time_t>(meta.build_timestamp);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::cout << "built: " << buf << "\n";
  }
  std::cout << "seed_extractor: " << seeds.extractor_name() << "\n";
  std::cout << "seed_k: " << get("k") << "\n";
  std::cout << "seed_window: " << get("window") << "\n";
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
  std::cout << "linearization: " << loaded.linearization_coords.size() << " nodes\n";
  std::cout << "1d_coordinates: " << (loaded.node_1d_coords.empty() ? "no" : "yes") << "\n";
  std::cout << "component_ids: " << (loaded.component_ids.empty() ? "no" : "yes") << "\n";

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
