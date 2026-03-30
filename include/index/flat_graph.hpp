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
 * - mmap-able: the layout can serve as a binary file format
 *
 * Memory per node: ~20 bytes (vs ~160 bytes for AlnNode)
 *
 * Related:
 *  - graph_store.hpp (abstract interface, FlatGraph plugs in via FlatGraphStore)
 *  - graph_store.hpp (abstract interface, FlatGraph plugs in)
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

  // Build from pre-assembled raw arrays (for deserialization)
  static FlatGraph fromRawArrays(
      std::uint32_t node_count, std::uint32_t path_count,
      std::vector<char> seq_data, std::vector<std::uint32_t> seq_offset,
      std::vector<std::uint32_t> seq_len,
      std::vector<char> name_data, std::vector<std::uint32_t> name_offset,
      std::vector<std::uint16_t> name_len,
      std::vector<std::uint8_t> is_reverse,
      std::vector<std::uint32_t> edge_target, std::vector<std::uint32_t> out_edge_offset,
      std::vector<std::uint32_t> step_data, std::vector<std::uint32_t> path_step_offset,
      std::vector<std::uint32_t> path_name_offset, std::vector<std::uint16_t> path_name_len,
      std::vector<std::uint64_t> path_length);

  // --- Node accessors ---

  std::uint32_t nodeCount() const { return node_count_; }

  std::string_view seq(std::uint32_t node_id) const {
    return {seq_data_.data() + seq_offset_[node_id], seq_len_[node_id]};
  }

  std::size_t seqLen(std::uint32_t node_id) const {
    return seq_len_[node_id];
  }

  std::string_view name(std::uint32_t node_id) const {
    return {name_data_.data() + name_offset_[node_id], name_len_[node_id]};
  }

  bool isReverse(std::uint32_t node_id) const {
    return is_reverse_[node_id] != 0;
  }

  // Label = original_id + direction suffix (e.g., "42+", "42-")
  std::string label(std::uint32_t node_id) const {
    auto n = name(node_id);
    return std::string(n) + (isReverse(node_id) ? "-" : "+");
  }

  // --- Edge accessors (CSR) ---

  std::uint32_t edgeCount() const { return static_cast<std::uint32_t>(edge_target_.size()); }

  // Outgoing neighbors of node_id
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

  std::uint64_t pathLength(std::uint32_t path_id) const {
    return path_length_[path_id];
  }

  void setPathLength(std::uint32_t path_id, std::uint64_t len) {
    path_length_[path_id] = len;
  }

  // Steps for path_id: array of node_ids
  const std::uint32_t* pathStepsBegin(std::uint32_t path_id) const {
    return step_data_.data() + path_step_offset_[path_id];
  }
  const std::uint32_t* pathStepsEnd(std::uint32_t path_id) const {
    return step_data_.data() + path_step_offset_[path_id + 1];
  }
  std::size_t pathStepCount(std::uint32_t path_id) const {
    return path_step_offset_[path_id + 1] - path_step_offset_[path_id];
  }

  // --- Raw arena access (for serialization) ---

  const std::vector<char>& seqData() const { return seq_data_; }
  const std::vector<std::uint32_t>& seqOffsets() const { return seq_offset_; }
  const std::vector<std::uint32_t>& seqLens() const { return seq_len_; }
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

  // Sequence arena
  std::vector<char> seq_data_;
  std::vector<std::uint32_t> seq_offset_;   // [node_count]
  std::vector<std::uint32_t> seq_len_;      // [node_count]

  // Name arena (shared by nodes and paths)
  std::vector<char> name_data_;
  std::vector<std::uint32_t> name_offset_;  // [node_count] for nodes
  std::vector<std::uint16_t> name_len_;     // [node_count] for nodes

  // Node metadata
  std::vector<std::uint8_t> is_reverse_;    // [node_count]

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
