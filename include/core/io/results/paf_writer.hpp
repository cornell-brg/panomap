// SPDX-License-Identifier: MIT
// Writer for PAF-formatted mapping results.

#pragma once

#include <fstream>
#include <string>

#include "core/index/flat_graph.hpp"
#include "core/io/results/result_writer.hpp"

namespace panomap::io {

class PafWriter : public ResultWriter {
public:
  PafWriter(const std::string& path, const index::FlatGraph& graph);
  ~PafWriter() override;

  void write(const mapping::ReadMapResult& result, const std::string& read_id,
             std::size_t read_length) override;

private:
  std::ofstream out_;
  const index::FlatGraph& graph_;
};

}  // namespace panomap::io
