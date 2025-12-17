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
    std::string query_sequence;
    std::string query_quality;

    char strand{'+'};  // '+' or '-'.

    std::string target_path;  // path/name in the reference graph (for PAF column 6).
    std::string graph_path;   // graph traversal like ">n1>n2>n3" (for GAF column 6).
    std::uint64_t target_length{0};
    std::uint64_t target_start{0};
    std::uint64_t target_end{0};

    std::uint64_t matches{0};
    std::uint64_t alignment_block_length{0};
    int mapq{0};

    struct Edit {
        std::uint32_t from_length{0};
        std::uint32_t to_length{0};
        std::string sequence;
    };

    struct Mapping {
        std::uint64_t node_id{0};
        std::uint32_t offset{0};
        bool is_reverse{false};
        std::vector<Edit> edits;
    };

    std::vector<Mapping> mappings;

    // Optional tags already formatted as GAF-style strings (e.g., "tp:A:P").
    std::vector<std::string> optional_fields;
};

}  // namespace piru::io
