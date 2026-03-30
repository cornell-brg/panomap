// SPDX-License-Identifier: MIT
// Factory helpers to select a result writer based on format/path.

#pragma once

#include <memory>
#include <string>

#include "index/flat_graph.hpp"
#include "io/results/result_writer.hpp"

namespace piru::io {

// Create a writer based on path extension.
ResultWriterPtr make_result_writer(const std::string& path,
                                    const piru::index::FlatGraph& graph);

// Create a writer with explicit format override.
ResultWriterPtr make_result_writer(const std::string& path,
                                    const std::string& format,
                                    const piru::index::FlatGraph& graph);

}  // namespace piru::io
