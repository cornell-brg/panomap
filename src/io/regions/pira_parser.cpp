// SPDX-License-Identifier: MIT

#include "io/regions/pira_parser.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include "util/logging.hpp"

namespace piru::io {

std::unordered_set<std::size_t> parse_pira(const std::string& filepath) {
  std::ifstream ifs(filepath);
  if (!ifs.is_open()) {
    throw std::runtime_error("parse_pira: cannot open '" + filepath + "'");
  }

  std::unordered_set<std::size_t> roi;
  std::string line;
  bool header_seen = false;

  while (std::getline(ifs, line)) {
    if (line.empty()) continue;

    // Header line
    if (line[0] == '#') {
      if (line.find("#pira") == 0) {
        header_seen = true;
        if (line.find("node-space=gfa") != std::string::npos) {
          throw std::runtime_error(
              "parse_pira: file uses GFA node IDs (--original-ids). "
              "Re-run `piru annotate` without --original-ids to get AlnGraph IDs.");
        }
        if (line.find("node-space=alngraph") == std::string::npos) {
          LOG_WARN("parse_pira: unrecognized node-space in header, assuming alngraph");
        }
      }
      continue;
    }

    if (!header_seen) {
      throw std::runtime_error("parse_pira: missing #pira header in '" + filepath + "'");
    }

    // Data line: "path:start-end id1,id2,id3,..."
    auto space_pos = line.find(' ');
    if (space_pos == std::string::npos || space_pos + 1 >= line.size()) continue;

    std::string ids_str = line.substr(space_pos + 1);
    std::size_t pos = 0;
    while (pos < ids_str.size()) {
      auto comma = ids_str.find(',', pos);
      if (comma == std::string::npos) comma = ids_str.size();
      std::string token = ids_str.substr(pos, comma - pos);
      if (!token.empty()) {
        roi.insert(std::stoull(token));
      }
      pos = comma + 1;
    }
  }

  LOG_INFO("parse_pira: loaded " + std::to_string(roi.size()) + " ROI node IDs from " + filepath);
  return roi;
}

}  // namespace piru::io
