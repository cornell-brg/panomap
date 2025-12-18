#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "io/results/result.hpp"
#include "io/results/result_writer_factory.hpp"
#include "mapping/chain_result_converter.hpp"
#include "mapping/seed_clusterer.hpp"

TEST_CASE("GAF writer writes basic record") {
    const auto tmp_path =
        std::filesystem::temp_directory_path() / "piru_test_output.gaf";

    auto writer = piru::io::make_result_writer(tmp_path.string());
    REQUIRE(writer != nullptr);

    piru::io::AlignmentResult r;
    r.query_name = "read1";
    r.query_length = 100;
    r.query_start = 10;
    r.query_end = 90;
    r.strand = '+';
    r.target_path = "chr1";
    r.target_length = 1000;
    r.target_start = 100;
    r.target_end = 180;
    r.matches = 70;
    r.alignment_block_length = 80;
    r.mapq = 60;
    r.optional_fields = {"tp:A:P"};

    REQUIRE(writer->write(r));
    writer.reset();

    std::ifstream in(tmp_path);
    REQUIRE(in.good());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    REQUIRE(lines.size() == 2);  // header + data
    CHECK(lines[0][0] == '#');   // header starts with #
    CHECK(lines[1] ==
          "read1\t100\t10\t90\t+\tchr1\t1000\t100\t180\t70\t80\t60\ttp:A:P");
}

TEST_CASE("Result writer factory produces json writer") {
    piru::io::AlignmentResult r;
    r.query_name = "q";
    r.query_length = 10;
    r.query_sequence = "ACGTACGTAA";
    r.target_path = "t";
    r.mappings.push_back({1, 0, false, {{4, 4, ""}}});

    auto json_writer = piru::io::make_result_writer("output.json");
    REQUIRE(json_writer != nullptr);
    CHECK(json_writer->write(r));
}

TEST_CASE("Result writer factory produces gam writer") {
    piru::io::AlignmentResult r;
    r.query_name = "q";
    r.query_length = 10;
    r.query_sequence = "ACGTACGTAA";
    r.target_path = "t";
    r.mappings.push_back({1, 0, false, {{4, 4, ""}}});

    auto gam_writer = piru::io::make_result_writer("output.gam");
    REQUIRE(gam_writer != nullptr);
    CHECK(gam_writer->write(r));
}

// =============================================================================
// DEV015: PAF Writer Tests
// =============================================================================

TEST_CASE("PAF writer writes basic record") {
    const auto tmp_path =
        std::filesystem::temp_directory_path() / "piru_test_output.paf";

    auto writer = piru::io::make_result_writer(tmp_path.string());
    REQUIRE(writer != nullptr);

    piru::io::AlignmentResult r;
    r.query_name = "read1";
    r.query_length = 100;
    r.query_start = 10;
    r.query_end = 90;
    r.strand = '+';
    r.target_path = "chr1";
    r.target_length = 1000;
    r.target_start = 100;
    r.target_end = 180;
    r.matches = 70;
    r.alignment_block_length = 80;
    r.mapq = 60;
    r.optional_fields = {"tp:A:P", "cs:i:150"};

    REQUIRE(writer->write(r));
    writer.reset();

    std::ifstream in(tmp_path);
    REQUIRE(in.good());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    REQUIRE(lines.size() == 2);  // header + data
    CHECK(lines[0][0] == '#');   // header starts with #
    // PAF: query_name, query_len, query_start, query_end, strand, target_name,
    //      target_len, target_start, target_end, matches, block_len, mapq, [tags]
    CHECK(lines[1] ==
          "read1\t100\t10\t90\t+\tchr1\t1000\t100\t180\t70\t80\t60\ttp:A:P\tcs:i:150");

    std::filesystem::remove(tmp_path);
}

TEST_CASE("PAF writer factory detects .paf extension") {
    auto writer = piru::io::make_result_writer("output.paf");
    REQUIRE(writer != nullptr);
}

// =============================================================================
// DEV015: ChainResultConverter Tests
// =============================================================================

namespace {

// Helper to create a simple test graph with paths
piru::index::AlnGraph makeTestGraph() {
    piru::index::AlnGraph graph;

    // Add nodes: node0 (s1), node1 (s2), node2 (s3)
    piru::index::AlnNode n0;
    n0.label = "s1";
    n0.original_id = "s1";
    n0.sequence = "ACGTACGT";
    n0.is_reverse = false;
    graph.addNode(n0);

    piru::index::AlnNode n1;
    n1.label = "s2";
    n1.original_id = "s2";
    n1.sequence = "GTACGTAC";
    n1.is_reverse = false;
    graph.addNode(n1);

    piru::index::AlnNode n2;
    n2.label = "s3";
    n2.original_id = "s3";
    n2.sequence = "TACGTACG";
    n2.is_reverse = false;
    graph.addNode(n2);

    // Add edges
    graph.addEdge({0, 1, 4});
    graph.addEdge({1, 2, 4});

    // Add a path: path1 -> s1 -> s2 -> s3
    piru::index::AlnPath path;
    path.name = "path1";
    path.steps = {{"s1", false}, {"s2", false}, {"s3", false}};
    path.overlaps = {4, 4};
    graph.addPath(path);

    return graph;
}

}  // namespace

TEST_CASE("ChainResultConverter converts single chain") {
    auto graph = makeTestGraph();
    piru::mapping::ChainResultConverter converter(graph);

    // Create a cluster summary with one chain
    piru::mapping::ClusterSummary summary;
    piru::mapping::ClusterGroup chain;
    chain.cluster_score = 100.0;

    // Add anchors spanning nodes 0, 1, 2
    piru::mapping::SeedAnchor a0;
    a0.target = {0, 2, 4};  // node_id=0, offset=2, length=4
    a0.read_pos = 10;
    a0.path_id = 0;
    a0.ref_coord = 2;
    chain.anchors.push_back(a0);

    piru::mapping::SeedAnchor a1;
    a1.target = {1, 0, 4};  // node_id=1, offset=0, length=4
    a1.read_pos = 14;
    a1.path_id = 0;
    a1.ref_coord = 6;
    chain.anchors.push_back(a1);

    piru::mapping::SeedAnchor a2;
    a2.target = {2, 0, 4};  // node_id=2, offset=0, length=4
    a2.read_pos = 18;
    a2.path_id = 0;
    a2.ref_coord = 10;
    chain.anchors.push_back(a2);

    summary.clusters.push_back(chain);

    // Convert
    auto results = converter.convert(summary, "test_read", 100);

    REQUIRE(results.size() == 1);
    CHECK(results[0].query_name == "test_read");
    CHECK(results[0].query_length == 100);
    CHECK(results[0].query_start == 10);
    CHECK(results[0].query_end == 22);  // 18 + 4 = 22
    CHECK(results[0].target_path == "path1");
    CHECK(results[0].target_start == 2);
    CHECK(results[0].target_end == 14);  // 10 + 4 = 14
    CHECK(results[0].strand == '+');
    CHECK(results[0].mapq == 60);  // No secondary, so high MAPQ
}

TEST_CASE("ChainResultConverter builds GAF path string") {
    auto graph = makeTestGraph();
    piru::mapping::ChainResultConverter converter(graph);

    piru::mapping::ClusterSummary summary;
    piru::mapping::ClusterGroup chain;
    chain.cluster_score = 50.0;

    // Anchors on nodes 0, 1, 2
    piru::mapping::SeedAnchor a0;
    a0.target = {0, 0, 4};
    a0.read_pos = 0;
    a0.path_id = 0;
    a0.ref_coord = 0;
    chain.anchors.push_back(a0);

    piru::mapping::SeedAnchor a1;
    a1.target = {1, 0, 4};
    a1.read_pos = 4;
    a1.path_id = 0;
    a1.ref_coord = 4;
    chain.anchors.push_back(a1);

    piru::mapping::SeedAnchor a2;
    a2.target = {2, 0, 4};
    a2.read_pos = 8;
    a2.path_id = 0;
    a2.ref_coord = 8;
    chain.anchors.push_back(a2);

    summary.clusters.push_back(chain);

    auto results = converter.convert(summary, "read1", 50);

    REQUIRE(results.size() == 1);
    // Path should be ">s1>s2>s3" (forward direction, original IDs)
    CHECK(results[0].graph_path == ">s1>s2>s3");
}

TEST_CASE("ChainResultConverter handles primary and secondary chains") {
    auto graph = makeTestGraph();
    piru::mapping::ChainResultConverter converter(graph);

    piru::mapping::ClusterSummary summary;

    // Primary chain (higher score)
    piru::mapping::ClusterGroup primary;
    primary.cluster_score = 100.0;
    piru::mapping::SeedAnchor a0;
    a0.target = {0, 0, 4};
    a0.read_pos = 0;
    a0.path_id = 0;
    a0.ref_coord = 0;
    primary.anchors.push_back(a0);
    summary.clusters.push_back(primary);

    // Secondary chain (lower score)
    piru::mapping::ClusterGroup secondary;
    secondary.cluster_score = 50.0;
    piru::mapping::SeedAnchor a1;
    a1.target = {1, 0, 4};
    a1.read_pos = 0;
    a1.path_id = 0;
    a1.ref_coord = 4;
    secondary.anchors.push_back(a1);
    summary.clusters.push_back(secondary);

    auto results = converter.convert(summary, "read1", 50);

    REQUIRE(results.size() == 2);

    // Primary should have higher MAPQ
    CHECK(results[0].mapq > 0);
    // Secondary should have MAPQ = 0
    CHECK(results[1].mapq == 0);

    // Check tp:A:P and tp:A:S tags
    bool primary_has_tag = false;
    bool secondary_has_tag = false;
    for (const auto& tag : results[0].optional_fields) {
        if (tag == "tp:A:P") primary_has_tag = true;
    }
    for (const auto& tag : results[1].optional_fields) {
        if (tag == "tp:A:S") secondary_has_tag = true;
    }
    CHECK(primary_has_tag);
    CHECK(secondary_has_tag);
}

TEST_CASE("ChainResultConverter respects primary_only config") {
    auto graph = makeTestGraph();
    piru::mapping::ChainResultConfig config;
    config.primary_only = true;
    piru::mapping::ChainResultConverter converter(graph, config);

    piru::mapping::ClusterSummary summary;

    // Add two chains
    piru::mapping::ClusterGroup c1;
    c1.cluster_score = 100.0;
    piru::mapping::SeedAnchor a0;
    a0.target = {0, 0, 4};
    a0.read_pos = 0;
    a0.path_id = 0;
    a0.ref_coord = 0;
    c1.anchors.push_back(a0);
    summary.clusters.push_back(c1);

    piru::mapping::ClusterGroup c2;
    c2.cluster_score = 50.0;
    piru::mapping::SeedAnchor a1;
    a1.target = {1, 0, 4};
    a1.read_pos = 0;
    a1.path_id = 0;
    a1.ref_coord = 4;
    c2.anchors.push_back(a1);
    summary.clusters.push_back(c2);

    auto results = converter.convert(summary, "read1", 50);

    // Should only get primary
    CHECK(results.size() == 1);
}

TEST_CASE("ChainResultConverter handles empty clusters") {
    auto graph = makeTestGraph();
    piru::mapping::ChainResultConverter converter(graph);

    piru::mapping::ClusterSummary summary;
    // No clusters

    auto results = converter.convert(summary, "read1", 50);

    CHECK(results.empty());
}
