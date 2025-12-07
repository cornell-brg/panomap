// SPDX-License-Identifier: MIT
// Pseudo-linearization utilities (SCC detection, chaining, coordinates).

#pragma once

#include <cstddef>
#include <vector>

#include "index/aln_graph.hpp"

namespace piru::index {

struct SccResult {
    std::vector<std::size_t> component;  // component id per node
    std::vector<bool> in_scc;            // true if node is in nontrivial SCC or self-loop
    std::size_t components{0};           // number of components
};

// Union-find data structures used across chaining steps.
struct UnionFind {
    std::vector<std::size_t> parent;
    std::vector<std::size_t> rank;

    explicit UnionFind(std::size_t n) : parent(n), rank(n, 0) {
        for (std::size_t i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    std::size_t find(std::size_t x);
    void unite(std::size_t x, std::size_t y);
};

struct TipFoldingResult {
    UnionFind uf;                        // Union-find structure with tips merged.
    std::vector<bool> ignorable_tip;     // True if node is part of a tip (skip in superbubble detection).

    explicit TipFoldingResult(std::size_t n) : uf(n), ignorable_tip(n, false) {}
};

TipFoldingResult chainTips(const AlnGraph& graph, const SccResult& scc);

void chainCycles(const AlnGraph& graph, TipFoldingResult& tip_result);

struct SuperbubbleResult {
    UnionFind uf;                    // Union-find structure with superbubbles merged.
    std::vector<bool> ignorable_tip; // Inherited from tip folding.

    SuperbubbleResult(UnionFind uf, std::vector<bool> ignorable_tip)
        : uf(std::move(uf)), ignorable_tip(std::move(ignorable_tip)) {}
};

SuperbubbleResult chainSuperbubbles(const AlnGraph& graph, const SccResult& scc,
                                     const TipFoldingResult& tip_result);

SccResult computeScc(const AlnGraph& graph);

}  // namespace piru::index
