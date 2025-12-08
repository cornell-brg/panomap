// SPDX-License-Identifier: MIT
#include <doctest/doctest.h>

#include "version.hpp"

#include <cstdio>
#include <filesystem>
#include <memory>

#include "index/aln_graph.hpp"
#include "index/graph_store.hpp"
#include "io/index/serialization.hpp"

namespace {
// Helper function to get a temporary file path
std::string temp_file_path(const std::string& name) {
    return std::filesystem::temp_directory_path() / name;
}
}

TEST_CASE("GraphStore serialization round-trip") {
    // 1. Create an in-memory AdjListGraphStore with some sample data.
    auto source_graph = std::make_unique<piru::index::AlnGraph>();
    piru::index::AlnNode n0, n1, n2;
    n0.label = "node0";
    n0.sequence = "GATTACA";
    n0.chain_id = 1;
    n0.linear_position = 100;
    source_graph->addNode(std::move(n0));

    n1.label = "node1";
    n1.sequence = "CAT";
    n1.chain_id = 1;
    n1.linear_position = 107;
    source_graph->addNode(std::move(n1));
    
    n2.label = "node2";
    n2.sequence = "TACC";
    n2.chain_id = 2;
    n2.linear_position = 50;
    source_graph->addNode(std::move(n2));
    
    source_graph->addEdge({0, 1, 3});
    source_graph->addEdge({1, 2, 1});

    piru::index::AdjListGraphStore source_store(std::move(*source_graph));

    // 2. Create some sample IndexMetadata.
    piru::io::index::IndexMetadata source_metadata;
    source_metadata.piru_version_major = 1;
    source_metadata.piru_version_minor = 2;
    source_metadata.piru_version_patch = 3;
    source_metadata.graph_flavor = 1; // DBG
    source_metadata.graph_k = 63;
    source_metadata.pore_k = 9;
    source_metadata.model_name = "test_model";
    source_metadata.fuzzy_quantizer = "rh2";
    source_metadata.align_quantizer = "int16";
    source_metadata.source_path = "test.gfa";

    const std::string test_path = temp_file_path("test.graph");

    // 3. Use write_graph to serialize the store.
    piru::io::index::write_graph(test_path, source_store, source_metadata);

    // 4. Use read_graph to deserialize it back.
    auto [loaded_store, loaded_metadata] = piru::io::index::read_graph(test_path);

    // 5. Assert that the deserialized store and metadata are identical.
    REQUIRE(loaded_store != nullptr);
    const auto& loaded_graph = loaded_store->graph();
    
    // Check metadata
    CHECK(loaded_metadata.piru_version_major == piru::Version::kMajor);
    CHECK(loaded_metadata.piru_version_minor == piru::Version::kMinor);
    CHECK(loaded_metadata.piru_version_patch == piru::Version::kPatch);
    CHECK(loaded_metadata.graph_flavor == 1);
    CHECK(loaded_metadata.graph_k == 63);
    CHECK(loaded_metadata.pore_k == 9);
    CHECK(loaded_metadata.model_name == "test_model");
    CHECK(loaded_metadata.fuzzy_quantizer == "rh2");
    CHECK(loaded_metadata.align_quantizer == "int16");
    CHECK(loaded_metadata.source_path == "test.gfa");

    // Check graph structure
    REQUIRE(loaded_graph.nodeCount() == 3);
    REQUIRE(loaded_graph.edgeCount() == 2);

    // Check node properties
    const auto& ln0 = loaded_graph.node(0);
    CHECK(ln0.label == "node0");
    CHECK(ln0.sequence == "GATTACA");
    CHECK(ln0.chain_id.value() == 1);
    CHECK(ln0.linear_position.value() == 100);

    const auto& ln1 = loaded_graph.node(1);
    CHECK(ln1.label == "node1");
    CHECK(ln1.sequence == "CAT");
    CHECK(ln1.chain_id.value() == 1);
    CHECK(ln1.linear_position.value() == 107);
    
    const auto& ln2 = loaded_graph.node(2);
    CHECK(ln2.label == "node2");
    CHECK(ln2.sequence == "TACC");
    CHECK(ln2.chain_id.value() == 2);
    CHECK(ln2.linear_position.value() == 50);

    // Check edges
    const auto& loaded_edges = loaded_graph.edges();
    bool edge0_found = false;
    bool edge1_found = false;
    for(const auto& edge : loaded_edges) {
        if(edge.from == 0 && edge.to == 1 && edge.overlap_bases == 3) {
            edge0_found = true;
        } else if (edge.from == 1 && edge.to == 2 && edge.overlap_bases == 1) {
            edge1_found = true;
        }
    }
    CHECK(edge0_found);
    CHECK(edge1_found);
    
    // Clean up the temporary file
    std::remove(test_path.c_str());
}
