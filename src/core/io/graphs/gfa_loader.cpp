#include "core/io/graphs/gfa_loader.hpp"

#include <fstream>
#include <sstream>
#include <vector>

#include "core/util/logging.hpp"

namespace piru::io {

namespace {

std::vector<std::string> split_tab(const std::string& line) {
  std::vector<std::string> parts;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    parts.push_back(field);
  }
  return parts;
}

std::vector<std::string> split_commas(const std::string& s) {
  std::vector<std::string> parts;
  std::stringstream ss(s);
  std::string field;
  while (std::getline(ss, field, ',')) {
    parts.push_back(field);
  }
  return parts;
}

bool parse_orientation(const std::string& token) { return !token.empty() && token[0] == '-'; }

std::optional<std::size_t> parse_overlap_bases(const std::string& cigar) {
  if (cigar.empty() || cigar == "*") return std::nullopt;
  std::size_t total = 0;
  std::size_t num = 0;
  for (char c : cigar) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      num = num * 10 + static_cast<std::size_t>(c - '0');
      continue;
    }
    // For overlap length, count operations that consume both sequences.
    if (c == 'M' || c == '=' || c == 'X') {
      total += num;
    }
    num = 0;
  }
  return total;
}

}  // namespace

GfaLoader::GfaLoader(std::string path) : path_(std::move(path)) {}

bool GfaLoader::load(ImportedGraph& graph) {
  graph.clear();
  std::ifstream in(path_);
  if (!in) {
    LOG_ERROR("Failed to open GFA file: " + path_);
    return false;
  }

  std::string line;
  std::size_t line_no = 0;
  bool ok = true;
  while (std::getline(in, line)) {
    ++line_no;
    if (line.empty() || line[0] == '#') continue;

    if (line.rfind("S\t", 0) == 0) {
      auto fields = split_tab(line);
      if (fields.size() < 3) {
        LOG_ERROR("Malformed GFA segment at line " + std::to_string(line_no));
        ok = false;
        continue;
      }
      graph.add_node({fields[1], fields[2]});
      continue;
    }

    if (line.rfind("L\t", 0) == 0) {
      auto fields = split_tab(line);
      if (fields.size() < 6) {
        LOG_ERROR("Malformed GFA link at line " + std::to_string(line_no));
        ok = false;
        continue;
      }
      ImportedGraphEdge edge;
      edge.from = fields[1];
      edge.from_reverse = parse_orientation(fields[2]);
      edge.to = fields[3];
      edge.to_reverse = parse_orientation(fields[4]);
      // Parse overlap CIGAR to numeric value. Assuming match-only overlaps (e.g., "0M").
      const std::string& overlap_field = fields[5];
      edge.overlap_bases = parse_overlap_bases(overlap_field);
      edge.overlap = edge.overlap_bases ? std::to_string(*edge.overlap_bases) : "0";
      graph.add_edge(std::move(edge));
      continue;
    }

    if (line.rfind("P\t", 0) == 0) {
      auto fields = split_tab(line);
      if (fields.size() < 4) {
        LOG_ERROR("Malformed GFA path at line " + std::to_string(line_no));
        ok = false;
        continue;
      }
      ImportedPath path;
      path.name = fields[1];

      const auto segments = split_commas(fields[2]);
      for (const auto& token : segments) {
        if (token.size() < 2) {
          LOG_ERROR("Malformed segment token in path '" + path.name + "' at line " +
                    std::to_string(line_no));
          ok = false;
          continue;
        }
        const char orient = token.back();
        if (orient != '+' && orient != '-') {
          LOG_ERROR("Missing orientation for segment '" + token + "' in path '" + path.name +
                    "' at line " + std::to_string(line_no));
          ok = false;
          continue;
        }
        ImportedPathStep step;
        step.segment_id = token.substr(0, token.size() - 1);
        step.is_reverse = (orient == '-');
        path.steps.push_back(std::move(step));
      }

      if (fields[3] != "*") {
        path.overlaps = split_commas(fields[3]);
        if (!path.steps.empty() && path.overlaps.size() + 1 != path.steps.size()) {
          LOG_WARN("Path '" + path.name + "' overlaps count mismatch; keeping raw overlaps");
        }
      }

      graph.add_path(std::move(path));
      continue;
    }

    // GFA 1.1 W (walk) lines. Format:
    //   W  SampleId  HapIndex  SeqId  SeqStart  SeqEnd  Walk
    // Walk uses '>' (forward) / '<' (reverse) prefix per segment, no commas.
    // Emit as a path named "Sample#Hap#Seq" (PanSN-spec) for piru's path index.
    if (line.rfind("W\t", 0) == 0) {
      auto fields = split_tab(line);
      if (fields.size() < 7) {
        LOG_ERROR("Malformed GFA walk at line " + std::to_string(line_no));
        ok = false;
        continue;
      }
      ImportedPath path;
      path.name = fields[1] + "#" + fields[2] + "#" + fields[3];

      const std::string& walk = fields[6];
      std::size_t i = 0;
      bool walk_ok = true;
      while (i < walk.size()) {
        char marker = walk[i];
        if (marker != '>' && marker != '<') {
          LOG_ERROR("Invalid walk char '" + std::string(1, marker) + "' in W line at " +
                    std::to_string(line_no));
          ok = false;
          walk_ok = false;
          break;
        }
        ++i;
        std::size_t start = i;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') ++i;
        if (i == start) {
          LOG_ERROR("Empty segment id in W walk at line " + std::to_string(line_no));
          ok = false;
          walk_ok = false;
          break;
        }
        ImportedPathStep step;
        step.segment_id = walk.substr(start, i - start);
        step.is_reverse = (marker == '<');
        path.steps.push_back(std::move(step));
      }
      if (walk_ok) graph.add_path(std::move(path));
      continue;
    }
  }

  return ok;
}

std::string GfaLoader::get_format_name() const { return "gfa"; }

}  // namespace piru::io
