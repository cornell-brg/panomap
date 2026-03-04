// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace piru::io {

struct BedRecord {
    std::string path_name;  // full GFA path name, e.g. SGDref#1#chrI
    std::int64_t start;     // 0-based inclusive
    std::int64_t end;       // 0-based exclusive
};

/// Parse a BED file. Each line: path_name<TAB>start<TAB>end [optional columns ignored].
/// Skips empty lines and lines starting with '#'.
std::vector<BedRecord> parse_bed(const std::string& filepath);

}  // namespace piru::io
