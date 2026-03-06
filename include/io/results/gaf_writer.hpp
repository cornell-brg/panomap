// SPDX-License-Identifier: MIT
// Writer for GAF-formatted alignment results.

#pragma once

#include <fstream>
#include <string>

#include "io/results/result_writer.hpp"

namespace piru::io {

class GafWriter : public ResultWriter {
public:
  explicit GafWriter(const std::string& path);
  ~GafWriter() override;

  bool write(const AlignmentResult& result) override;

private:
  std::ofstream out_;
};

using GafWriterPtr = std::unique_ptr<GafWriter>;

}  // namespace piru::io
