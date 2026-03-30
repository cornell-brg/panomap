// SPDX-License-Identifier: MIT

#include "index/flat_graph.hpp"
#include "index/aln_graph.hpp"

namespace piru::index {

FlatGraph FlatGraph::fromAlnGraph(const AlnGraph& aln) {
  FlatGraph fg;
  fg.node_count_ = static_cast<std::uint32_t>(aln.nodeCount());
  fg.path_count_ = static_cast<std::uint32_t>(aln.pathCount());

  // --- Nodes: sequences + names + metadata ---

  // Pre-compute total sequence and name lengths for reservation
  std::size_t total_seq = 0;
  std::size_t total_name = 0;
  for (std::size_t i = 0; i < aln.nodeCount(); ++i) {
    total_seq += aln.node(i).sequence.size();
    total_name += aln.node(i).original_id.size();
  }
  for (std::size_t i = 0; i < aln.pathCount(); ++i) {
    total_name += aln.paths()[i].name.size();
  }

  fg.seq_data_.reserve(total_seq);
  fg.name_data_.reserve(total_name);
  fg.seq_offset_.resize(fg.node_count_);
  fg.seq_len_.resize(fg.node_count_);
  fg.name_offset_.resize(fg.node_count_);
  fg.name_len_.resize(fg.node_count_);
  fg.is_reverse_.resize(fg.node_count_);

  for (std::uint32_t i = 0; i < fg.node_count_; ++i) {
    const auto& node = aln.node(i);

    // Sequence
    fg.seq_offset_[i] = static_cast<std::uint32_t>(fg.seq_data_.size());
    fg.seq_len_[i] = static_cast<std::uint32_t>(node.sequence.size());
    fg.seq_data_.insert(fg.seq_data_.end(), node.sequence.begin(), node.sequence.end());

    // Name (use original_id which is what serialization uses)
    fg.name_offset_[i] = static_cast<std::uint32_t>(fg.name_data_.size());
    fg.name_len_[i] = static_cast<std::uint16_t>(node.original_id.size());
    fg.name_data_.insert(fg.name_data_.end(), node.original_id.begin(), node.original_id.end());

    // Metadata
    fg.is_reverse_[i] = node.is_reverse ? 1 : 0;
  }

  // --- Edges (CSR) ---

  fg.out_edge_offset_.resize(fg.node_count_ + 1);
  std::size_t total_edges = 0;
  for (std::uint32_t i = 0; i < fg.node_count_; ++i) {
    fg.out_edge_offset_[i] = static_cast<std::uint32_t>(total_edges);
    total_edges += aln.outgoing(i).size();
  }
  fg.out_edge_offset_[fg.node_count_] = static_cast<std::uint32_t>(total_edges);

  fg.edge_target_.resize(total_edges);
  for (std::uint32_t i = 0; i < fg.node_count_; ++i) {
    const auto& out = aln.outgoing(i);
    std::uint32_t offset = fg.out_edge_offset_[i];
    for (std::size_t j = 0; j < out.size(); ++j) {
      fg.edge_target_[offset + j] = static_cast<std::uint32_t>(out[j]);
    }
  }

  // --- Paths (CSR) ---

  fg.path_step_offset_.resize(fg.path_count_ + 1);
  fg.path_name_offset_.resize(fg.path_count_);
  fg.path_name_len_.resize(fg.path_count_);
  fg.path_length_.resize(fg.path_count_);

  std::size_t total_steps = 0;
  for (std::uint32_t i = 0; i < fg.path_count_; ++i) {
    fg.path_step_offset_[i] = static_cast<std::uint32_t>(total_steps);
    total_steps += aln.paths()[i].steps.size();
  }
  fg.path_step_offset_[fg.path_count_] = static_cast<std::uint32_t>(total_steps);

  fg.step_data_.resize(total_steps);
  for (std::uint32_t i = 0; i < fg.path_count_; ++i) {
    const auto& path = aln.paths()[i];

    // Path name (into shared name arena)
    fg.path_name_offset_[i] = static_cast<std::uint32_t>(fg.name_data_.size());
    fg.path_name_len_[i] = static_cast<std::uint16_t>(path.name.size());
    fg.name_data_.insert(fg.name_data_.end(), path.name.begin(), path.name.end());

    // Path length
    fg.path_length_[i] = path.length;

    // Steps (node_ids)
    std::uint32_t offset = fg.path_step_offset_[i];
    for (std::size_t j = 0; j < path.steps.size(); ++j) {
      fg.step_data_[offset + j] = static_cast<std::uint32_t>(std::stoull(path.steps[j].node_id));
    }
  }

  return fg;
}

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
