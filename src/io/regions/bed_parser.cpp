// SPDX-License-Identifier: MIT

#include "io/regions/bed_parser.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "util/logging.hpp"

namespace piru::io {

std::vector<BedRecord> parse_bed(const std::string& filepath) {
  std::ifstream ifs(filepath);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open BED file: " + filepath);
  }

  std::vector<BedRecord> records;
  std::string line;
  int line_num = 0;

  while (std::getline(ifs, line)) {
    ++line_num;

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream ss(line);
    std::string path_name;
    std::int64_t start, end;

    if (!(ss >> path_name >> start >> end)) {
      LOG_WARN("BED line " + std::to_string(line_num) + ": skipping malformed line");
      continue;
    }

    if (start < 0 || end < 0 || end <= start) {
      LOG_WARN("BED line " + std::to_string(line_num) + ": invalid interval [" +
               std::to_string(start) + ", " + std::to_string(end) + ")");
      continue;
    }

    records.push_back({std::move(path_name), start, end});
  }

  LOG_INFO("Parsed " + std::to_string(records.size()) + " BED records from " + filepath);
  return records;
}

}  // namespace piru::io
