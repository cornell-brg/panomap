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
// Different implementations handle different linearization strategies:
// - PathWalkExpander: Uses linearization_coords (1 hit → N anchors)
// - SuperbubbleExpander: Uses GraphStore chain_id (1 hit → 1 anchor)
class AnchorExpander {
public:
    virtual ~AnchorExpander() = default;

    // Expand seed hits to anchors.
    // Returns vector of anchors (may be larger than input if nodes appear on multiple paths).
    virtual std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const = 0;

    // Backend name for logging and debugging.
    virtual std::string name() const = 0;
};

// Path-walk expansion: uses linearization coordinates to expand seed hits.
//
// For each seed hit:
// - Lookup linearization coordinates for the node
// - Emit one anchor per coordinate occurrence (handles multiple paths/cycles)
// - Add node_offset to ref_coord to get precise position
//
// Multiplicity: 1 hit → N anchors (N = number of paths the node appears on)
// Seeds with no linearization coordinates (empty vector) are skipped.
class PathWalkExpander : public AnchorExpander {
public:
    // Construct expander with linearization coordinates from a Linearizer.
    // coords[node_id] = vector of LinearCoordinate for that node
    explicit PathWalkExpander(const std::vector<std::vector<index::LinearCoordinate>>& coords);

    std::vector<Anchor> expand(const std::vector<SeedHitRecord>& hits) const override;

    std::string name() const override { return "path-walk"; }

private:
    const std::vector<std::vector<index::LinearCoordinate>>& coords_;
};

// Superbubble expansion: uses GraphStore chain_id and linear_position.
//
// For each seed hit:
// - Lookup chain_id and linear_position from GraphStore
// - Create anchor with path_id = chain_id, ref_coord = linear_position + offset
//
// Multiplicity: 1 hit → 1 anchor (trivial mapping)
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
