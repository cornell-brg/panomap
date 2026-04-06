// SPDX-License-Identifier: MIT
// Factory helpers to select a result writer based on format/path.

#pragma once

#include <memory>
#include <string>

#include "index/flat_graph.hpp"
#include "io/results/result_writer.hpp"

namespace piru::io {

// Create a writer based on path extension (.paf -> PAF, anything else -> GAF).
ResultWriterPtr make_result_writer(const std::string& path, const piru::index::FlatGraph& graph,
                                   bool primary_only = true);

// Create a GAF writer to stdout.
ResultWriterPtr make_result_writer_stdout(const piru::index::FlatGraph& graph,
                                          bool primary_only = true);

}  // namespace piru::io
