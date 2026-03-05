// SPDX-License-Identifier: MIT
// Index pipeline tests
// Note: pseudo_linearize and superbubble_linearizer tests removed in DEV040

#include <doctest/doctest.h>

#include "index/aln_graph.hpp"

using namespace piru;

// Placeholder - add index pipeline integration tests as needed
TEST_CASE("AlnGraph basic construction") {
    index::AlnGraph g;
    g.addNode({});
    g.addNode({});
    g.addEdge({0, 1, 0});

    CHECK(g.nodeCount() == 2);
    CHECK(g.edgeCount() == 1);
}
