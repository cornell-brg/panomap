// SPDX-License-Identifier: MIT
// Expands seed hits to anchors using linearization coordinates.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "mapping/anchor.hpp"
#include "mapping/seed_clusterer.hpp"

namespace piru::mapping {

// Abstract interface for expanding seed hits (graph space) to anchors (linear space).
//
// This is the key transformation stage that converts graph-based coordinates (node_id, offset)
// into linear reference coordinates (path_id, ref_coord) to enable efficient clustering/chaining.
//
// Different implementations handle different linearization strategies:
// - SuperbubbleExpander: Trivial 1:1 mapping using GraphStore chain_id (superbubble pipeline)
// - PathWalkExpander: 1:N mapping via path occurrence coordinates (path-walk pipeline)
//
// When to use:
// - Superbubble: Simple variation graphs, local coordinates within bubbles, fast O(n) clustering
// - Path-walk: Complex graphs with cycles, haplotype-aware mapping, colinear chaining
class AnchorExpander {
public:
    virtual ~AnchorExpander() = default;

    // Expand seed hits to anchors.
    // Returns vector of anchors (may be larger than input if nodes appear on multiple paths).
    virtual std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const = 0;

    // Backend name for logging and debugging.
    virtual std::string name() const = 0;
};

using AnchorExpanderPtr = std::unique_ptr<AnchorExpander>;

// Path-walk expansion: uses linearization coordinates to expand seed hits.
//
// For each seed hit:
// - Lookup linearization coordinates for the node
// - Emit one anchor per coordinate occurrence (handles multiple paths/cycles)
// - Add node_offset to ref_coord to get precise position
//
// Multiplicity: 1 hit -> N anchors (N = number of paths the node appears on)
// Seeds with no linearization coordinates (empty vector) are skipped.
class PathWalkExpander : public AnchorExpander {
public:
    // Construct expander with linearization coordinates and path lengths.
    // coords[node_id] = vector of LinearCoordinate for that node
    // path_lengths[path_id] = length of that path (for bounds checking)
    PathWalkExpander(const std::vector<std::vector<index::LinearCoordinate>>& coords,
                     const std::vector<std::size_t>& path_lengths);

    std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const override;

    std::string name() const override { return "path-walk"; }

private:
    const std::vector<std::vector<index::LinearCoordinate>>& coords_;
    const std::vector<std::size_t>& path_lengths_;
};

// Superbubble expansion: uses GraphStore chain_id and linear_position.
//
// For each seed hit:
// - Lookup chain_id and linear_position from GraphStore
// - Create anchor with path_id = chain_id, ref_coord = linear_position + offset
//
// Multiplicity: 1 hit -> 1 anchor (trivial mapping)
// Seeds without chain_id (unmapped nodes) are skipped.
class SuperbubbleExpander : public AnchorExpander {
public:
    // Construct expander with GraphStore containing chain_id and linear_position per node.
    explicit SuperbubbleExpander(const index::GraphStore* graph_store);

    std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const override;

    std::string name() const override { return "superbubble"; }

private:
    const index::GraphStore* graph_store_;
};

}  // namespace piru::mapping
