// SPDX-License-Identifier: MIT
// Minimal alignment result representation for output writers.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace piru::io {

struct AlignmentResult {
    std::string query_name;
    std::uint64_t query_length{0};
    std::uint64_t query_start{0};
    std::uint64_t query_end{0};

    char strand{'+'};  // '+' or '-'.

    std::string target_path;  // path/name in the reference graph.
    std::uint64_t target_length{0};
    std::uint64_t target_start{0};
    std::uint64_t target_end{0};

    std::uint64_t matches{0};
    std::uint64_t alignment_block_length{0};
    int mapq{0};

    // Optional tags already formatted as GAF-style strings (e.g., "tp:A:P").
    std::vector<std::string> optional_fields;
};

}  // namespace piru::io
