// SPDX-License-Identifier: MIT
// Writer for PAF-formatted mapping results.

#pragma once

#include <fstream>
#include <string>

#include "core/io/results/result_writer.hpp"

namespace piru::io {

class PafWriter : public ResultWriter {
public:
  explicit PafWriter(const std::string& path);
  ~PafWriter() override;

  void write(const mapping::ReadMapResult& result, const std::string& read_id,
             std::size_t read_length) override;

private:
  std::ofstream out_;
};

}  // namespace piru::io
