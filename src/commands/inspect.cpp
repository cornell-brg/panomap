/**
 * inspect.cpp
 *
 * CLI handler for `piru inspect`. Loads a .pirx index and prints
 * metadata, graph stats, and seed store stats.
 *
 * SPDX-License-Identifier: MIT
 */

#include "commands/inspect.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "io/index/serialization.hpp"
#include "util/logging.hpp"

int handle_inspect(const std::vector<std::string>& args) {
  piru::cli::Parsed parsed;
  piru::cli::ParseConfig config;
  config.usage = "Usage: piru inspect <index.pirx>";
  config.positional_help = {"<index.pirx>    Index file to inspect"};
  config.options = {
      {'h', "help", false, "Show help"},
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

  const auto& meta = loaded.metadata;
  const auto& fg = loaded.graph->flat();
  const auto& seeds = *loaded.seeds;

  const auto& params = seeds.params();
  auto get = [&](const std::string& key) -> std::string {
    auto it = params.find(key);
    return it != params.end() ? it->second : "";
  };

  std::cout << "index: " << index_path << "\n";
  std::cout << "piru_version: " << (meta.version.empty() ? "unknown" : meta.version) << "\n";
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

  return 0;
}
