// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>

namespace piru::io {

/// Parse a .pira annotation file (produced by `piru annotate`).
/// Returns a flat set of all ROI node IDs (AlnGraph space) across all intervals.
/// Errors if the file uses GFA node space (--original-ids).
std::unordered_set<std::size_t> parse_pira(const std::string& filepath);

}  // namespace piru::io
