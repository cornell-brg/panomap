// SPDX-License-Identifier: MIT
// Writer for PAF-formatted alignment results.

#pragma once

#include <fstream>
#include <string>

#include "io/results/result_writer.hpp"

namespace piru::io {

class PafWriter : public ResultWriter {
public:
  explicit PafWriter(const std::string& path);
  ~PafWriter() override;

  bool write(const AlignmentResult& result) override;

private:
  std::ofstream out_;
};

using PafWriterPtr = std::unique_ptr<PafWriter>;

}  // namespace piru::io
