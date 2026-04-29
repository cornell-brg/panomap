// SPDX-License-Identifier: MIT

#include "core/index/flat_graph.hpp"

#include <cstring>

namespace piru::index {

FlatGraph FlatGraph::fromRawArrays(
    std::uint32_t node_count, std::uint32_t path_count, std::vector<char> seq_data,
    std::vector<std::uint32_t> seq_offset, std::vector<std::uint32_t> seq_len,
    std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
    std::vector<std::uint16_t> name_len, std::vector<std::uint8_t> is_reverse,
    std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
    std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
    std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
    std::vector<std::uint64_t> path_length) {
  FlatGraph fg;
  fg.node_count_ = node_count;
  fg.path_count_ = path_count;
  fg.seq_len_ = std::move(seq_len);

  // Two construction modes:
  //  1. BUILD path (e.g., test fixtures, GFA loaders): seq_data is non-empty,
  //     contains ASCII bases. Pack into 2-bit + N mask, populate
  //     seq_base_offset_ + seq_packed_ + seq_n_mask_ + total_bases_.
  //  2. LOAD path (deserialization from .pirx): seq_data is EMPTY -- the
  //     packed arena is build-time-only and isn't serialized to disk.
  //     Skip allocations entirely; only the metadata (names, lens, edges,
  //     paths) is needed at map time. base2bit / isN must NOT be called
  //     after this construction (they're only used by indexers, not the
  //     mapper). On a 100M-node pangenome this saves ~800MB at map time.
  //
  // **DO NOT** restore unconditional allocation here without auditing
  // base2bit / isN call sites -- there's a uint32 overflow in the packed
  // offsets above 4GB cumulative bases (fix dev-108 c9ba503).
  if (seq_data.empty()) {
    // Load-only mode: leave seq_base_offset_, seq_packed_, seq_n_mask_ empty.
    // total_bases_ stays 0; setSeqLens() may overwrite seq_len_ later but
    // does not need to reconstruct the packed arena.
  } else {
    fg.seq_base_offset_.resize(node_count);

    // Compute total bases and base offsets
    fg.total_bases_ = 0;
    for (std::uint32_t i = 0; i < node_count; ++i) {
      fg.seq_base_offset_[i] = fg.total_bases_;
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
      std::uint64_t base_off = fg.seq_base_offset_[i];
      std::uint32_t len = fg.seq_len_[i];

      for (std::uint32_t j = 0; j < len; ++j) {
        char c = src[j];
        std::uint64_t abs_pos = base_off + j;

        // Encode base
        std::uint8_t val = encode2bit(c);
        std::uint64_t byte_idx = abs_pos >> 2;
        std::uint32_t bit_shift = (abs_pos & 3) << 1;
        fg.seq_packed_[byte_idx] |= (val << bit_shift);

        // Mark N positions
        if (c == 'N' || c == 'n') {
          fg.seq_n_mask_[abs_pos >> 3] |= (1 << (abs_pos & 7));
        }
      }
    }
  }
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

FlatGraph FlatGraph::fromPackedArrays(
    std::uint32_t node_count, std::uint32_t path_count, std::size_t total_bases,
    std::vector<std::uint8_t> seq_packed, std::vector<std::uint8_t> seq_n_mask,
    std::vector<std::uint64_t> seq_base_offset, std::vector<std::uint32_t> seq_len,
    std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
    std::vector<std::uint16_t> name_len, std::vector<std::uint8_t> is_reverse,
    std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
    std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
    std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
    std::vector<std::uint64_t> path_length) {
  FlatGraph fg;
  fg.node_count_ = node_count;
  fg.path_count_ = path_count;
  fg.total_bases_ = total_bases;
  fg.seq_packed_ = std::move(seq_packed);
  fg.seq_n_mask_ = std::move(seq_n_mask);
  fg.seq_base_offset_ = std::move(seq_base_offset);
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
