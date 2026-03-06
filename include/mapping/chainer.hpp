// SPDX-License-Identifier: MIT
// Interfaces for chaining backends.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/seed_store.hpp"

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

// Minimal hit record used for clustering/chaining.
struct NodeAnchor {
  index::SeedEntry target;    // node_id + offset in graph
  std::size_t read_pos{0};    // seed position in the read
  std::uint64_t hash{0};      // seed hash (for debugging/uniqueness)
  std::size_t span{0};        // coverage length on query (from Seed.length, may be merged)
  std::size_t frequency{0};   // occurrences of this hash in the index
  mutable double score{0.0};  // Computed during clustering (mutable for legacy compatibility)
};

// Anchor candidate produced by clustering/chaining.
struct ChainedAnchor {
  index::SeedEntry target;  // node_id + offset in graph
  std::size_t read_pos{0};  // position in read
  double score{0.0};        // backend-specific score
  std::size_t chain_id{0};  // which chain this anchor belongs to

  // Optional: linear coordinates (for path-walk pipeline debugging)
  std::size_t path_id{0};     // reference path ID
  std::int64_t ref_coord{0};  // linear position on reference path
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
