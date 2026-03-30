// SPDX-License-Identifier: MIT

#include "index/flat_graph.hpp"

namespace piru::index {

FlatGraph FlatGraph::fromRawArrays(
    std::uint32_t node_count, std::uint32_t path_count,
    std::vector<char> seq_data, std::vector<std::uint32_t> seq_offset,
    std::vector<std::uint32_t> seq_len,
    std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
    std::vector<std::uint16_t> name_len,
    std::vector<std::uint8_t> is_reverse,
    std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
    std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
    std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
    std::vector<std::uint64_t> path_length) {
  FlatGraph fg;
  fg.node_count_ = node_count;
  fg.path_count_ = path_count;
  fg.seq_data_ = std::move(seq_data);
  fg.seq_offset_ = std::move(seq_offset);
  fg.seq_len_ = std::move(seq_len);
  fg.name_data_ = std::move(name_data);
  fg.name_offset_ = std::move(name_offset);
  fg.name_len_ = std::move(name_len);
  fg.is_reverse_ = std::move(is_reverse);
  fg.edge_target_ = std::move(edge_target);
  fg.out_edge_offset_ = std::move(out_edge_offset);
  fg.step_data_ = std::move(step_data);
  fg.path_step_offset_ = std::move(path_step_offset);
  fg.path_name_offset_ = std::move(path_name_offset);
  fg.path_name_len_ = std::move(path_name_len);
  fg.path_length_ = std::move(path_length);
  return fg;
}

}  // namespace piru::index
