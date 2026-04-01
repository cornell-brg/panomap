// SPDX-License-Identifier: MIT
// Writer for GAF-formatted mapping results.

#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include "index/flat_graph.hpp"
#include "io/results/result_writer.hpp"

namespace piru::io {

struct GafWriterConfig {
  bool primary_only{false};
  std::size_t max_secondary{9};
  double min_secondary_ratio{0.7};
};

class GafWriter : public ResultWriter {
public:
  GafWriter(const std::string& path, const index::FlatGraph& graph, GafWriterConfig config = {});
  ~GafWriter() override;

  void write(const mapping::ReadMapResult& result, const std::string& read_id,
             std::size_t read_length) override;

private:
  // Build GAF-style node walk string (">n1>n2>n3")
  std::string buildPathString(const std::vector<mapping::ChainedAnchor>& anchors) const;

  // Get path name and length from graph
  std::string getPathName(std::size_t path_id) const;
  std::size_t getPathLength(std::size_t path_id) const;

  std::ofstream out_;
  const index::FlatGraph& graph_;
  GafWriterConfig config_;

  mutable std::vector<std::size_t> path_lengths_;
  mutable bool path_lengths_computed_{false};
};

}  // namespace piru::io
