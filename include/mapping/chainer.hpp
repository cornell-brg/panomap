// SPDX-License-Identifier: MIT
// Interfaces for chaining backends.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cli/parse.hpp"

namespace piru::mapping {

// Minimal hit record: graph-space seed hit (16 bytes).
// Produced by SeedLookup, consumed by chainer and seed merger.
struct NodeAnchor {
  std::uint32_t node_id{0};   // Graph node ID
  std::uint32_t offset{0};    // Offset within node
  std::uint32_t read_pos{0};  // Seed position in the read
  std::uint16_t span{0};      // Coverage length on query
  std::uint16_t length{0};    // Coverage length on reference (may differ after merge)
};

// Anchor candidate produced by chaining (24B).
// score/chain_id removed (dead), path_id lifted to Chain.
struct ChainedAnchor {
  std::uint32_t node_id{0};   // Graph node ID
  std::uint32_t offset{0};    // Offset within node
  std::uint32_t read_pos{0};  // Position in read
  std::uint16_t length{0};    // Coverage length
  std::int64_t ref_coord{0};  // Linear position on reference path
};

// Coordinate space for chain ref_coords.
enum class CoordSpace {
  kPath,      // ref_coord in path-linear space (PathChainer, GraphChainer)
  kCanonical, // ref_coord in 1D canonical space (SortChainer)
};

// A single chain: scored group of anchors.
struct Chain {
  double score{0.0};
  std::size_t path_id{0};                     // Reference path ID (kPath only)
  CoordSpace coord_space{CoordSpace::kPath};  // Coordinate system for ref_coords
  std::vector<ChainedAnchor> anchors;
};

struct ChainResult {
  std::vector<Chain> chains;
  std::size_t expanded_anchor_count{0};  // total anchors before chaining
  std::vector<bool> used_inputs;         // used_inputs[i] = true if input hits[i] participated in a chain
};

// Abstract interface for chaining backends.
// Backends receive NodeAnchors (graph-space) and decide internally how to
// transform and chain them (e.g. expand to linear space, chain in graph space).
class Chainer {
public:
  virtual ~Chainer() = default;

  virtual ChainResult chain(const std::vector<NodeAnchor>& hits) const = 0;
  virtual std::string name() const = 0;
};

using ChainerPtr = std::unique_ptr<Chainer>;

}  // namespace piru::mapping
