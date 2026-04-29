#pragma once

/**
 * flat_graph.hpp
 *
 * FlatGraph: memory-efficient graph representation using flat arenas and
 * CSR adjacency. Replaces AlnGraph's heap-allocated strings and nested vectors.
 *
 * Design principles (from FlatGFA / Adrian Sampson, Cornell):
 * - All data in contiguous arrays (no inner heap allocations)
 * - u32 indices instead of pointers (2x savings on 64-bit)
 * - CSR for adjacency lists (no vector-of-vector)
 * - Arena for variable-length data (sequences, names)
 * - 2-bit packed DNA sequences (4x compression vs ASCII)
 *
 * Related:
 *  - graph_store.hpp (abstract interface, FlatGraph plugs in via FlatGraphStore)
 *
 * SPDX-License-Identifier: MIT
 */

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace piru::index {

class FlatGraph {
public:
  FlatGraph() = default;

  // Build from pre-assembled raw arrays (ASCII sequences -- packed internally)
  static FlatGraph fromRawArrays(
      std::uint32_t node_count, std::uint32_t path_count, std::vector<char> seq_data,
      std::vector<std::uint32_t> seq_offset, std::vector<std::uint32_t> seq_len,
      std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
      std::vector<std::uint16_t> name_len, std::vector<std::uint8_t> is_reverse,
      std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
      std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
      std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
      std::vector<std::uint64_t> path_length);

  // Build from pre-packed 2-bit arrays (no ASCII intermediate)
  static FlatGraph fromPackedArrays(
      std::uint32_t node_count, std::uint32_t path_count, std::size_t total_bases,
      std::vector<std::uint8_t> seq_packed, std::vector<std::uint8_t> seq_n_mask,
      std::vector<std::uint64_t> seq_base_offset, std::vector<std::uint32_t> seq_len,
      std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
      std::vector<std::uint16_t> name_len, std::vector<std::uint8_t> is_reverse,
      std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
      std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
      std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
      std::vector<std::uint64_t> path_length);

  // --- 2-bit encoding helpers ---

  // Encode: A=0, C=1, G=2, T=3. N and others map to 0 (check n_mask).
  static constexpr std::uint8_t encode2bit(char c) {
    // Lookup table approach for speed
    switch (c) {
      case 'C':
      case 'c':
        return 1;
      case 'G':
      case 'g':
        return 2;
      case 'T':
      case 't':
        return 3;
      default:
        return 0;  // A, N, and anything else
    }
  }

  static constexpr char decode2bit(std::uint8_t v) {
    constexpr char table[] = {'A', 'C', 'G', 'T'};
    return table[v & 3];
  }

  // --- Node accessors ---

  std::uint32_t nodeCount() const { return node_count_; }

  std::size_t seqLen(std::uint32_t node_id) const { return seq_len_[node_id]; }

  // Get 2-bit encoded base at position within node (0-3, check isN for ambiguity).
  // abs_pos / byte_idx must be 64-bit: total_bases on hg38 (fwd+rev) is ~6.2 GB,
  // overflowing uint32. (Bug fixed dev-108 -- silent corruption past 4 GB cumulative.)
  std::uint8_t base2bit(std::uint32_t node_id, std::uint32_t pos) const {
    std::uint64_t abs_pos = seq_base_offset_[node_id] + pos;
    std::uint64_t byte_idx = abs_pos >> 2;         // abs_pos / 4
    std::uint32_t bit_shift = (abs_pos & 3) << 1;  // (abs_pos % 4) * 2
    return (seq_packed_[byte_idx] >> bit_shift) & 3;
  }

  // Check if position is an N base
  bool isN(std::uint32_t node_id, std::uint32_t pos) const {
    std::uint64_t abs_pos = seq_base_offset_[node_id] + pos;
    return (seq_n_mask_[abs_pos >> 3] >> (abs_pos & 7)) & 1;
  }

  // Decode full node sequence to ASCII string (for serialization/export)
  std::string seqDecoded(std::uint32_t node_id) const {
    std::size_t len = seq_len_[node_id];
    std::string result(len, '\0');
    for (std::size_t i = 0; i < len; ++i) {
      if (isN(node_id, static_cast<std::uint32_t>(i))) {
        result[i] = 'N';
      } else {
        result[i] = decode2bit(base2bit(node_id, static_cast<std::uint32_t>(i)));
      }
    }
    return result;
  }

  std::string_view name(std::uint32_t node_id) const {
    return {name_data_.data() + name_offset_[node_id], name_len_[node_id]};
  }

  bool isReverse(std::uint32_t node_id) const { return is_reverse_[node_id] != 0; }

  // Label = original_id + direction suffix (e.g., "42+", "42-")
  std::string label(std::uint32_t node_id) const {
    auto n = name(node_id);
    return std::string(n) + (isReverse(node_id) ? "-" : "+");
  }

  // --- Edge accessors (CSR) ---

  std::uint32_t edgeCount() const { return static_cast<std::uint32_t>(edge_target_.size()); }

  const std::uint32_t* outBegin(std::uint32_t node_id) const {
    return edge_target_.data() + out_edge_offset_[node_id];
  }
  const std::uint32_t* outEnd(std::uint32_t node_id) const {
    return edge_target_.data() + out_edge_offset_[node_id + 1];
  }
  std::size_t outDegree(std::uint32_t node_id) const {
    return out_edge_offset_[node_id + 1] - out_edge_offset_[node_id];
  }

  // --- Path accessors (CSR) ---

  std::uint32_t pathCount() const { return path_count_; }

  std::string_view pathName(std::uint32_t path_id) const {
    return {name_data_.data() + path_name_offset_[path_id], path_name_len_[path_id]};
  }

  std::uint64_t pathLength(std::uint32_t path_id) const { return path_length_[path_id]; }

  void setPathLength(std::uint32_t path_id, std::uint64_t len) { path_length_[path_id] = len; }

  const std::uint32_t* pathStepsBegin(std::uint32_t path_id) const {
    return step_data_.data() + path_step_offset_[path_id];
  }
  const std::uint32_t* pathStepsEnd(std::uint32_t path_id) const {
    return step_data_.data() + path_step_offset_[path_id + 1];
  }
  std::size_t pathStepCount(std::uint32_t path_id) const {
    return path_step_offset_[path_id + 1] - path_step_offset_[path_id];
  }

  // --- Raw access (for serialization) ---

  const std::vector<std::uint8_t>& seqPacked() const { return seq_packed_; }
  const std::vector<std::uint8_t>& seqNMask() const { return seq_n_mask_; }
  const std::vector<std::uint64_t>& seqBaseOffsets() const { return seq_base_offset_; }
  const std::vector<std::uint32_t>& seqLens() const { return seq_len_; }
  void setSeqLens(std::vector<std::uint32_t> lens) { seq_len_ = std::move(lens); }
  std::size_t totalBases() const { return total_bases_; }
  const std::vector<char>& nameData() const { return name_data_; }
  const std::vector<std::uint32_t>& nameOffsets() const { return name_offset_; }
  const std::vector<std::uint16_t>& nameLens() const { return name_len_; }
  const std::vector<std::uint8_t>& isReverseVec() const { return is_reverse_; }
  const std::vector<std::uint32_t>& edgeTargets() const { return edge_target_; }
  const std::vector<std::uint32_t>& outEdgeOffsets() const { return out_edge_offset_; }
  const std::vector<std::uint32_t>& stepData() const { return step_data_; }
  const std::vector<std::uint32_t>& pathStepOffsets() const { return path_step_offset_; }
  const std::vector<std::uint16_t>& pathNameLens() const { return path_name_len_; }
  const std::vector<std::uint32_t>& pathNameOffsets() const { return path_name_offset_; }
  const std::vector<std::uint64_t>& pathLengths() const { return path_length_; }

private:
  std::uint32_t node_count_{0};
  std::uint32_t path_count_{0};
  std::size_t total_bases_{0};

  // 2-bit packed sequence arena (4 bases per byte)
  std::vector<std::uint8_t> seq_packed_;        // ceil(total_bases / 4) bytes
  std::vector<std::uint8_t> seq_n_mask_;        // 1 bit per base (ceil(total_bases / 8) bytes)
  std::vector<std::uint64_t> seq_base_offset_;  // [node_count] offset in bases
  std::vector<std::uint32_t> seq_len_;          // [node_count] length in bases

  // Name arena (shared by nodes and paths)
  std::vector<char> name_data_;
  std::vector<std::uint32_t> name_offset_;  // [node_count] for nodes
  std::vector<std::uint16_t> name_len_;     // [node_count] for nodes

  // Node metadata
  std::vector<std::uint8_t> is_reverse_;  // [node_count]

  // Edges (CSR: out_edge_offset[node_count+1], edge_target[edge_count])
  std::vector<std::uint32_t> edge_target_;
  std::vector<std::uint32_t> out_edge_offset_;  // [node_count + 1]

  // Paths (CSR: path_step_offset[path_count+1], step_data[total_steps])
  std::vector<std::uint32_t> step_data_;
  std::vector<std::uint32_t> path_step_offset_;  // [path_count + 1]
  std::vector<std::uint32_t> path_name_offset_;  // [path_count] into name_data
  std::vector<std::uint16_t> path_name_len_;     // [path_count]
  std::vector<std::uint64_t> path_length_;       // [path_count] in base pairs
};

}  // namespace piru::index
