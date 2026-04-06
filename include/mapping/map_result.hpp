/**
 * map_result.hpp
 *
 * Mapping result types for the mapping pipeline. Mapping holds a single
 * chain result, ReadMapResult holds all mappings for a read.
 *
 * Related:
 *  - batch_mapper.cpp (produces ReadMapResult)
 *  - gaf_writer.cpp (consumes ReadMapResult)
 *  - chainer.hpp (ChainedAnchor, CoordSpace)
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <vector>

#include "mapping/chainer.hpp"

namespace piru::mapping {

/**
 * A single mapping (one chain).
 *
 * Populated by the chainer, then MAPQ computed in batch_mapper.
 */
struct Mapping {
  std::vector<ChainedAnchor> anchors;
  double chain_score{0.0};
  std::size_t path_id{0};                     // Reference path ID (kPath only)
  CoordSpace coord_space{CoordSpace::kPath};  // Coordinate system for ref_coords
  int mapq{0};                                // Mapping quality (computed post-chaining)
};

/**
 * All mappings for a single read.
 */
struct ReadMapResult {
  std::vector<Mapping> mappings;  // primary (index 0) + secondaries
  std::size_t expanded_anchor_count{0};

  /* Per-read timing */
  std::size_t chunks_processed{0};
  double processing_time_sec{0.0};

  bool mapped() const { return !mappings.empty(); }
  const Mapping* primary() const { return mappings.empty() ? nullptr : &mappings[0]; }
};

}  // namespace piru::mapping
