// SPDX-License-Identifier: MIT
// Writer for GAF-formatted mapping results.

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/index/flat_graph.hpp"
#include "core/io/results/result_writer.hpp"

namespace panomap::io {

struct GafWriterConfig {
  bool primary_only{true};
  std::size_t max_secondary{9};
  const std::vector<float>* node_1d_coords{nullptr};
  const std::vector<std::uint32_t>* component_ids{nullptr};
  // List only the nodes the anchors fall on, with path positions in columns 7 to 9,
  // instead of the full walk with walk positions.
  bool anchor_nodes{false};
};

class GafWriter : public ResultWriter {
public:
  // Write to file
  GafWriter(const std::string& path, const index::FlatGraph& graph, GafWriterConfig config = {});
  // Write to stream (e.g. stdout)
  GafWriter(const index::FlatGraph& graph, GafWriterConfig config = {});
  ~GafWriter() override;

  void write(const mapping::ReadMapResult& result, const std::string& read_id,
             std::size_t read_length) override;

private:
  std::string buildPathString(const std::vector<mapping::ChainedAnchor>& anchors) const;
  // Start of each step on its path, plus the path length at the end. Built on first use.
  const std::vector<std::uint64_t>& stepStarts(std::uint32_t path_id) const;

  std::ofstream file_;        // owned file (when writing to file)
  std::ostream* out_;         // points to file_ or external stream (e.g. cout)
  const index::FlatGraph& graph_;
  GafWriterConfig config_;

  mutable std::unordered_map<std::uint32_t, std::vector<std::uint64_t>> step_starts_;
};

}  // namespace panomap::io
