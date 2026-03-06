#include "io/results/gaf_writer.hpp"

#include <sstream>

#include "util/logging.hpp"

namespace piru::io {

namespace {

std::string to_gaf_line(const AlignmentResult& r) {
  std::stringstream ss;
  // GAF column 6 is graph path traversal (e.g., ">n1>n2>n3"), not path name
  const std::string& path_col = r.graph_path.empty() ? r.target_path : r.graph_path;
  ss << r.query_name << '\t' << r.query_length << '\t' << r.query_start << '\t' << r.query_end
     << '\t' << r.strand << '\t' << path_col << '\t' << r.target_length << '\t' << r.target_start
     << '\t' << r.target_end << '\t' << r.matches << '\t' << r.alignment_block_length << '\t'
     << r.mapq;
  // Add path name as pn:Z: tag if we have a graph path
  if (!r.graph_path.empty() && !r.target_path.empty()) {
    ss << '\t' << "pn:Z:" << r.target_path;
  }
  for (const auto& tag : r.optional_fields) {
    ss << '\t' << tag;
  }
  return ss.str();
}

}  // namespace

GafWriter::GafWriter(const std::string& path) : out_(path, std::ios::out | std::ios::trunc) {
  if (!out_) {
    LOG_ERROR("Failed to open GAF output: " + path);
  } else {
    // Write header comment for development convenience
    out_ << "#query_name\tquery_len\tquery_start\tquery_end\tstrand\t"
         << "graph_path\tpath_len\tpath_start\tpath_end\t"
         << "matches\tblock_len\tmapq\t[optional_tags: pn:Z:path_name ...]\n";
  }
}

GafWriter::~GafWriter() { out_.flush(); }

bool GafWriter::write(const AlignmentResult& result) {
  if (!out_) return false;
  out_ << to_gaf_line(result) << '\n';
  return static_cast<bool>(out_);
}

}  // namespace piru::io
