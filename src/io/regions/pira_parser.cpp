// SPDX-License-Identifier: MIT

#include "io/regions/pira_parser.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include "util/logging.hpp"

namespace piru::io {

namespace {

void parse_node_list(const std::string& ids_str, std::unordered_set<std::size_t>& out) {
  std::size_t pos = 0;
  while (pos < ids_str.size()) {
    auto comma = ids_str.find(',', pos);
    if (comma == std::string::npos) comma = ids_str.size();
    std::string token = ids_str.substr(pos, comma - pos);
    if (!token.empty()) {
      out.insert(std::stoull(token));
    }
    pos = comma + 1;
  }
}

}  // namespace

PiraFile parse_pira_v2(const std::string& filepath) {
  std::ifstream ifs(filepath);
  if (!ifs.is_open()) {
    throw std::runtime_error("parse_pira: cannot open '" + filepath + "'");
  }

  PiraFile result;
  std::string line;
  bool is_v2 = false;

  while (std::getline(ifs, line)) {
    if (line.empty()) continue;

    if (line[0] == '#') {
      if (line.find("#pira v2") != std::string::npos) is_v2 = true;
      continue;
    }

    PiraRegion region;

    if (is_v2) {
      // v2: label\t1d_start\t1d_end\tnodes
      auto tab1 = line.find('\t');
      if (tab1 == std::string::npos) continue;
      auto tab2 = line.find('\t', tab1 + 1);
      if (tab2 == std::string::npos) continue;
      auto tab3 = line.find('\t', tab2 + 1);
      if (tab3 == std::string::npos) continue;

      region.label = line.substr(0, tab1);
      std::string s_start = line.substr(tab1 + 1, tab2 - tab1 - 1);
      std::string s_end = line.substr(tab2 + 1, tab3 - tab2 - 1);

      if (s_start != "*" && s_end != "*") {
        region.coord_1d_start = std::stod(s_start);
        region.coord_1d_end = std::stod(s_end);
        region.has_1d = true;
      }

      parse_node_list(line.substr(tab3 + 1), region.nodes);
    } else {
      // v1: label node1,node2,...
      auto space_pos = line.find(' ');
      if (space_pos == std::string::npos) continue;
      region.label = line.substr(0, space_pos);
      parse_node_list(line.substr(space_pos + 1), region.nodes);
    }

    for (std::size_t nid : region.nodes) {
      result.all_nodes.insert(nid);
    }
    result.regions.push_back(std::move(region));
  }

  LOG_INFO("parse_pira: " + std::to_string(result.regions.size()) + " regions, " +
           std::to_string(result.all_nodes.size()) + " total nodes");
  return result;
}

std::unordered_set<std::size_t> parse_pira(const std::string& filepath) {
  auto pira = parse_pira_v2(filepath);
  return std::move(pira.all_nodes);
}

}  // namespace piru::io
