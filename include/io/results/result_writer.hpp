// SPDX-License-Identifier: MIT
// Interface for writing mapping results.

#pragma once

#include <memory>
#include <string>

#include "mapping/map_result.hpp"

namespace piru::io {

class ResultWriter {
public:
  virtual ~ResultWriter() = default;

  // Write results for a single read.
  virtual void write(const mapping::ReadMapResult& result, const std::string& read_id,
                     std::size_t read_length) = 0;
};

using ResultWriterPtr = std::unique_ptr<ResultWriter>;

}  // namespace piru::io
