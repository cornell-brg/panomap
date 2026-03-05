// SPDX-License-Identifier: MIT

#include "io/results/paf_writer.hpp"

#include <sstream>

#include "util/logging.hpp"

namespace piru::io {

namespace {

// PAF format: 12 mandatory columns + optional tags
// 1. Query name
// 2. Query length
// 3. Query start (0-based)
// 4. Query end
// 5. Strand (+/-)
// 6. Target name
// 7. Target length
// 8. Target start (0-based)
// 9. Target end
// 10. Number of matches
// 11. Alignment block length
// 12. Mapping quality
std::string to_paf_line(const AlignmentResult& r) {
    std::stringstream ss;
    ss << r.query_name << '\t' << r.query_length << '\t' << r.query_start << '\t' << r.query_end
       << '\t' << r.strand << '\t' << r.target_path << '\t' << r.target_length << '\t'
       << r.target_start << '\t' << r.target_end << '\t' << r.matches << '\t'
       << r.alignment_block_length << '\t' << r.mapq;

    // Optional tags
    for (const auto& tag : r.optional_fields) {
        ss << '\t' << tag;
    }

    return ss.str();
}

}  // namespace

PafWriter::PafWriter(const std::string& path) : out_(path, std::ios::out | std::ios::trunc) {
    if (!out_) {
        LOG_ERROR("Failed to open PAF output: " + path);
    } else {
        // Write header comment for development convenience
        out_ << "#query_name\tquery_len\tquery_start\tquery_end\tstrand\t"
             << "target_path\ttarget_len\ttarget_start\ttarget_end\t"
             << "matches\tblock_len\tmapq\t[optional_tags...]\n";
    }
}

PafWriter::~PafWriter() { out_.flush(); }

bool PafWriter::write(const AlignmentResult& result) {
    if (!out_) return false;
    out_ << to_paf_line(result) << '\n';
    return static_cast<bool>(out_);
}

}  // namespace piru::io
