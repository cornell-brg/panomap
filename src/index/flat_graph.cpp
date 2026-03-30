// SPDX-License-Identifier: MIT

#include "index/flat_graph.hpp"

#include <cstring>

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

  // Pack ASCII sequences into 2-bit encoding + N mask
  // seq_data contains ASCII chars, seq_offset/seq_len index into it
  fg.seq_len_ = std::move(seq_len);
  fg.seq_base_offset_.resize(node_count);

  // Compute total bases and base offsets
  fg.total_bases_ = 0;
  for (std::uint32_t i = 0; i < node_count; ++i) {
    fg.seq_base_offset_[i] = static_cast<std::uint32_t>(fg.total_bases_);
    fg.total_bases_ += fg.seq_len_[i];
  }

  // Allocate packed arrays
  std::size_t packed_bytes = (fg.total_bases_ + 3) / 4;  // ceil(total / 4)
  std::size_t mask_bytes = (fg.total_bases_ + 7) / 8;    // ceil(total / 8)
  fg.seq_packed_.resize(packed_bytes, 0);
  fg.seq_n_mask_.resize(mask_bytes, 0);

  // Pack each node's sequence
  for (std::uint32_t i = 0; i < node_count; ++i) {
    const char* src = seq_data.data() + seq_offset[i];
    std::uint32_t base_off = fg.seq_base_offset_[i];
    std::uint32_t len = fg.seq_len_[i];

    for (std::uint32_t j = 0; j < len; ++j) {
      char c = src[j];
      std::uint32_t abs_pos = base_off + j;

      // Encode base
      std::uint8_t val = encode2bit(c);
      std::uint32_t byte_idx = abs_pos >> 2;
      std::uint32_t bit_shift = (abs_pos & 3) << 1;
      fg.seq_packed_[byte_idx] |= (val << bit_shift);

      // Mark N positions
      if (c == 'N' || c == 'n') {
        fg.seq_n_mask_[abs_pos >> 3] |= (1 << (abs_pos & 7));
      }
    }
  }

  // Release the ASCII sequence data (no longer needed)
  // seq_data and seq_offset are consumed, not stored

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
