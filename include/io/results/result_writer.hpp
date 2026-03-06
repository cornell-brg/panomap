// SPDX-License-Identifier: MIT
// Interface for writing alignment results.

#pragma once

#include <memory>

#include "io/results/result.hpp"

namespace piru::io {

class ResultWriter {
public:
  virtual ~ResultWriter() = default;

  // Write a single result. Returns false on failure.
  virtual bool write(const AlignmentResult& result) = 0;
};

using ResultWriterPtr = std::unique_ptr<ResultWriter>;

}  // namespace piru::io
