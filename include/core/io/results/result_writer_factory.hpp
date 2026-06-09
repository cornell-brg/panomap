// SPDX-License-Identifier: MIT
// Factory helpers to select a result writer based on format/path.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "core/index/flat_graph.hpp"
#include "core/io/results/result_writer.hpp"

namespace panomap::io {

// Create a writer based on path extension (.paf -> PAF, anything else -> GAF).
ResultWriterPtr make_result_writer(const std::string& path, const panomap::index::FlatGraph& graph,
                                   bool primary_only = true,
                                   const std::vector<float>* node_1d_coords = nullptr,
                                   const std::vector<std::uint32_t>* component_ids = nullptr);

// Create a GAF writer to stdout.
ResultWriterPtr make_result_writer_stdout(const panomap::index::FlatGraph& graph,
                                          bool primary_only = true,
                                          const std::vector<float>* node_1d_coords = nullptr,
                                          const std::vector<std::uint32_t>* component_ids = nullptr);

}  // namespace panomap::io
