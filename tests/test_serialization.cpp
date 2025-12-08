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

TEST_CASE("SignalStore serialization round-trip (int16)") {
    // 1. Create a VectorSignalStore with int16 data.
    std::vector<piru::signal::AlignmentQuantizedSignal> signals;
    signals.emplace_back(piru::signal::AlignmentQuantizedSignal{
        .kind = piru::signal::AlignmentQuantizationKind::kInt16,
        .data = std::vector<int16_t>{100, 101, 102}
    });
    signals.emplace_back(piru::signal::AlignmentQuantizedSignal{
        .kind = piru::signal::AlignmentQuantizationKind::kInt16,
        .data = std::vector<int16_t>{-50, 0, 50, 100}
    });
    piru::index::VectorSignalStore source_store(std::move(signals));
    
    const std::string test_path = temp_file_path("test_int16.signals");

    // 2. Use write_signals to serialize the store.
    const float test_scale = 0.5f;
    const float test_offset = 50.0f;
    piru::io::index::write_signals(test_path, source_store, test_scale, test_offset);

    // 3. Use read_signals to deserialize it back.
    auto [loaded_store, loaded_metadata] = piru::io::index::read_signals(test_path);

    // 4. Assert that the deserialized store and metadata are identical.
    REQUIRE(loaded_store != nullptr);
    CHECK(loaded_metadata.quantization_bits == 16);
    CHECK(loaded_metadata.scale == test_scale);
    CHECK(loaded_metadata.offset == test_offset);
    REQUIRE(loaded_store->size() == 2);

    const auto* sig0 = loaded_store->get(0);
    REQUIRE(sig0 != nullptr);
    CHECK(sig0->kind == piru::signal::AlignmentQuantizationKind::kInt16);
    const auto& data0 = std::get<std::vector<int16_t>>(sig0->data);
    CHECK(data0.size() == 3);
    CHECK(data0[0] == 100);
    CHECK(data0[1] == 101);
    CHECK(data0[2] == 102);

    const auto* sig1 = loaded_store->get(1);
    REQUIRE(sig1 != nullptr);
    CHECK(sig1->kind == piru::signal::AlignmentQuantizationKind::kInt16);
    const auto& data1 = std::get<std::vector<int16_t>>(sig1->data);
    CHECK(data1.size() == 4);
    CHECK(data1[0] == -50);
    CHECK(data1[1] == 0);
    CHECK(data1[2] == 50);
    CHECK(data1[3] == 100);

    // Clean up
    std::remove(test_path.c_str());
}

TEST_CASE("SignalStore serialization round-trip (float32)") {
    // 1. Create a VectorSignalStore with float32 data.
    std::vector<piru::signal::AlignmentQuantizedSignal> signals;
    signals.emplace_back(piru::signal::AlignmentQuantizedSignal{
        .kind = piru::signal::AlignmentQuantizationKind::kFloat32,
        .data = std::vector<float>{100.5f, 101.25f, 102.0f}
    });
    piru::index::VectorSignalStore source_store(std::move(signals));

    const std::string test_path = temp_file_path("test_float32.signals");

    // 2. Use write_signals to serialize the store.
    piru::io::index::write_signals(test_path, source_store, 1.0f, 0.0f);

    // 3. Use read_signals to deserialize it back.
    auto [loaded_store, loaded_metadata] = piru::io::index::read_signals(test_path);

    // 4. Assert that the deserialized store is identical.
    REQUIRE(loaded_store != nullptr);
    CHECK(loaded_metadata.quantization_bits == 32);
    CHECK(loaded_metadata.scale == 1.0f);
    CHECK(loaded_metadata.offset == 0.0f);
    REQUIRE(loaded_store->size() == 1);

    const auto* sig0 = loaded_store->get(0);
    REQUIRE(sig0 != nullptr);
    CHECK(sig0->kind == piru::signal::AlignmentQuantizationKind::kFloat32);
    const auto& data0 = std::get<std::vector<float>>(sig0->data);
    CHECK(data0.size() == 3);
    CHECK(data0[0] == 100.5f);
    CHECK(data0[1] == 101.25f);
    CHECK(data0[2] == 102.0f);

    // Clean up
    std::remove(test_path.c_str());
}

TEST_CASE("SignalStore serialization round-trip (int8)") {
    // 1. Create a VectorSignalStore with int8 data.
    std::vector<piru::signal::AlignmentQuantizedSignal> signals;
    signals.emplace_back(piru::signal::AlignmentQuantizedSignal{
        .kind = piru::signal::AlignmentQuantizationKind::kInt8,
        .data = std::vector<int8_t>{-10, 0, 10, 20}
    });
    piru::index::VectorSignalStore source_store(std::move(signals));

    const std::string test_path = temp_file_path("test_int8.signals");

    // 2. Use write_signals to serialize the store.
    piru::io::index::write_signals(test_path, source_store, 1.0f, 0.0f);

    // 3. Use read_signals to deserialize it back.
    auto [loaded_store, loaded_metadata] = piru::io::index::read_signals(test_path);

    // 4. Assert that the deserialized store is identical.
    REQUIRE(loaded_store != nullptr);
    REQUIRE(loaded_store->size() == 1);
    CHECK(loaded_metadata.quantization_bits == 8);

    const auto* sig0 = loaded_store->get(0);
    REQUIRE(sig0 != nullptr);
    CHECK(sig0->kind == piru::signal::AlignmentQuantizationKind::kInt8);
    const auto& data0 = std::get<std::vector<int8_t>>(sig0->data);
    CHECK(data0.size() == 4);
    CHECK(data0[0] == -10);
    CHECK(data0[1] == 0);
    CHECK(data0[2] == 10);
    // Clean up
    std::remove(test_path.c_str());
}

TEST_CASE("SeedStore serialization round-trip") {
    // 1. Create a HashSeedStore with sample data.
    auto source_store = std::make_unique<piru::index::HashSeedStore>();
    source_store->insert(123, {1, 10});
    source_store->insert(456, {2, 20});
    source_store->insert(123, {3, 30});
    source_store->set_extractor_name("kmer");
    source_store->set_params({{"k", "10"}, {"stride", "5"}});
    source_store->set_max_hash_frequency(2);
    source_store->set_frequency_threshold(100);
    source_store->set_filter_fraction(0.1);

    const std::string test_path = temp_file_path("test.seeds");

    // 2. Use write_seeds to serialize it.
    piru::io::index::write_seeds(test_path, *source_store);

    // 3. Use read_seeds to deserialize it.
    auto loaded_store = piru::io::index::read_seeds(test_path);

    // 4. Assert that the stores are identical.
    REQUIRE(loaded_store != nullptr);
    CHECK(loaded_store->size() == 2);
    CHECK(loaded_store->extractor_name() == "kmer");
    const auto& params = loaded_store->params();
    REQUIRE(params.size() == 2);
    CHECK(params.at("k") == "10");
    CHECK(params.at("stride") == "5");
    CHECK(loaded_store->max_hash_frequency() == 2);
    CHECK(loaded_store->frequency_threshold() == 100);
    CHECK(loaded_store->filter_fraction() == 0.1);
    
    const auto* hits123 = loaded_store->lookup(123);
    REQUIRE(hits123 != nullptr);
    REQUIRE(hits123->size() == 2);
    CHECK((*hits123)[0].node_id == 1);
    CHECK((*hits123)[0].offset == 10);
    CHECK((*hits123)[1].node_id == 3);
    CHECK((*hits123)[1].offset == 30);
    
    const auto* hits456 = loaded_store->lookup(456);
    REQUIRE(hits456 != nullptr);
    REQUIRE(hits456->size() == 1);
    CHECK((*hits456)[0].node_id == 2);
    CHECK((*hits456)[0].offset == 20);

    const auto* hits999 = loaded_store->lookup(999);
    CHECK(hits999 == nullptr);

    // Clean up
    std::remove(test_path.c_str());
}
