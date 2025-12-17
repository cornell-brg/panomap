// SPDX-License-Identifier: MIT
// Factory helpers to select a result writer based on format/path.

#pragma once

#include <memory>
#include <string>

#include "io/results/result_writer.hpp"

namespace piru::io {

// Attempt to choose a writer based on path/extension.
// Returns nullptr if no suitable writer is available.
ResultWriterPtr make_result_writer(const std::string& path);

// Create a writer with explicit format override.
// Format should be one of: "paf", "gaf", "gam", "json"
// Returns nullptr if format is unsupported.
ResultWriterPtr make_result_writer(const std::string& path, const std::string& format);

}  // namespace piru::io
