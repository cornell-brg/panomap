#include "io/results/gaf_writer.hpp"

#include <sstream>

#include "util/logging.hpp"

namespace piru::io {

namespace {

std::string to_gaf_line(const AlignmentResult& r) {
    std::stringstream ss;
    ss << r.query_name << '\t' << r.query_length << '\t' << r.query_start << '\t'
       << r.query_end << '\t' << r.strand << '\t' << r.target_path << '\t'
       << r.target_length << '\t' << r.target_start << '\t' << r.target_end << '\t'
       << r.matches << '\t' << r.alignment_block_length << '\t' << r.mapq;
    for (const auto& tag : r.optional_fields) {
        ss << '\t' << tag;
    }
    return ss.str();
}

}  // namespace

GafWriter::GafWriter(const std::string& path) : out_(path, std::ios::out | std::ios::trunc) {
    if (!out_) {
        LOG_ERROR("Failed to open GAF output: " + path);
    }
}

GafWriter::~GafWriter() { out_.flush(); }

bool GafWriter::write(const AlignmentResult& result) {
    if (!out_) return false;
    out_ << to_gaf_line(result) << '\n';
    return static_cast<bool>(out_);
}

}  // namespace piru::io
