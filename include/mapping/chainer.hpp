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

// Compact anchor in linear path coordinate space (16 bytes).
// Produced by expanding NodeAnchors using linearization coordinates.
// One NodeAnchor may expand to multiple PathAnchors (one per path occurrence).
// Back-references to graph space are recovered via src_idx into the NodeAnchor
// input array -- keeps this struct small for cache efficiency in the DP loop.
struct PathAnchor {
  std::uint32_t ref_coord{0};  // Linear position along reference path
  std::uint32_t query_pos{0};  // Position in query/read
  std::uint16_t length{0};     // Coverage length (from seed span)
  std::uint16_t _pad{0};       // Reserved
  std::uint32_t src_idx{0};    // Index into NodeAnchor input array
};

// Minimal hit record: graph-space seed hit (16 bytes).
// Produced by SeedLookup, consumed by chainer and seed merger.
struct NodeAnchor {
  std::uint32_t node_id{0};   // Graph node ID
  std::uint32_t offset{0};    // Offset within node
  std::uint32_t read_pos{0};  // Seed position in the read
  std::uint16_t span{0};      // Coverage length on query
  std::uint16_t length{0};    // Coverage length on reference (may differ after merge)
};

// Anchor candidate produced by chaining.
struct ChainedAnchor {
  std::uint32_t node_id{0};   // Graph node ID
  std::uint32_t offset{0};    // Offset within node
  std::uint32_t read_pos{0};  // Position in read
  std::uint16_t length{0};    // Coverage length
  std::uint16_t _pad{0};      // Reserved
  double score{0.0};          // Backend-specific score
  std::size_t chain_id{0};    // Which chain this anchor belongs to
  std::size_t path_id{0};     // Reference path ID
  std::int64_t ref_coord{0};  // Linear position on reference path
};

// A single chain: scored group of anchors.
struct Chain {
  double score{0.0};
  std::vector<ChainedAnchor> anchors;
};

struct ChainResult {
  double score{0.0};
  std::vector<ChainedAnchor> anchors;    // flat list from best chain
  std::vector<Chain> chains;             // all extracted chains
  std::size_t expanded_anchor_count{0};  // total anchors before chaining
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
