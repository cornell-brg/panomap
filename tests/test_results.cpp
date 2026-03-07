#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "io/results/result.hpp"
#include "io/results/result_writer_factory.hpp"
#include "mapping/map_result.hpp"
#include "mapping/result_formatter.hpp"

TEST_CASE("GAF writer writes basic record") {
  const auto tmp_path = std::filesystem::temp_directory_path() / "piru_test_output.gaf";

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
  CHECK(lines[1] == "read1\t100\t10\t90\t+\tchr1\t1000\t100\t180\t70\t80\t60\ttp:A:P");
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
  const auto tmp_path = std::filesystem::temp_directory_path() / "piru_test_output.paf";

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
  CHECK(lines[1] == "read1\t100\t10\t90\t+\tchr1\t1000\t100\t180\t70\t80\t60\ttp:A:P\tcs:i:150");

  std::filesystem::remove(tmp_path);
}

TEST_CASE("PAF writer factory detects .paf extension") {
  auto writer = piru::io::make_result_writer("output.paf");
  REQUIRE(writer != nullptr);
}

// =============================================================================
// DEV015: ResultFormatter Tests
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

TEST_CASE("ResultFormatter formats single mapping") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatter formatter(graph);

  // Create a ReadMapResult with one mapping
  piru::mapping::ReadMapResult map_result;
  piru::mapping::Mapping mapping;
  mapping.chain_score = 100.0;

  // Add anchors spanning nodes 0, 1, 2
  piru::mapping::ChainedAnchor a0;
  a0.node_id = 0; a0.offset = 2; a0.length = 4;
  a0.read_pos = 10;
  a0.path_id = 0;
  a0.ref_coord = 2;
  mapping.anchors.push_back(a0);

  piru::mapping::ChainedAnchor a1;
  a1.node_id = 1; a1.offset = 0; a1.length = 4;
  a1.read_pos = 14;
  a1.path_id = 0;
  a1.ref_coord = 6;
  mapping.anchors.push_back(a1);

  piru::mapping::ChainedAnchor a2;
  a2.node_id = 2; a2.offset = 0; a2.length = 4;
  a2.read_pos = 18;
  a2.path_id = 0;
  a2.ref_coord = 10;
  mapping.anchors.push_back(a2);

  map_result.mappings.push_back(mapping);

  // Format
  auto results = formatter.format(map_result, "test_read", 100);

  REQUIRE(results.size() == 1);
  CHECK(results[0].query_name == "test_read");
  CHECK(results[0].query_length == 100);
  CHECK(results[0].query_start == 10);
  CHECK(results[0].query_end == 22);  // 18 + 4 = 22
  CHECK(results[0].target_path == "path1");
  CHECK(results[0].target_start == 2);
  CHECK(results[0].target_end == 14);  // 10 + 4 = 14
  CHECK(results[0].strand == '+');
  CHECK(results[0].mapq == 100);  // MAPQ equals raw chain score
}

TEST_CASE("ResultFormatter builds GAF path string") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatter formatter(graph);

  piru::mapping::ReadMapResult map_result;
  piru::mapping::Mapping mapping;
  mapping.chain_score = 50.0;

  // Anchors on nodes 0, 1, 2
  piru::mapping::ChainedAnchor a0;
  a0.node_id = 0; a0.offset = 0; a0.length = 4;
  a0.read_pos = 0;
  a0.path_id = 0;
  a0.ref_coord = 0;
  mapping.anchors.push_back(a0);

  piru::mapping::ChainedAnchor a1;
  a1.node_id = 1; a1.offset = 0; a1.length = 4;
  a1.read_pos = 4;
  a1.path_id = 0;
  a1.ref_coord = 4;
  mapping.anchors.push_back(a1);

  piru::mapping::ChainedAnchor a2;
  a2.node_id = 2; a2.offset = 0; a2.length = 4;
  a2.read_pos = 8;
  a2.path_id = 0;
  a2.ref_coord = 8;
  mapping.anchors.push_back(a2);

  map_result.mappings.push_back(mapping);

  auto results = formatter.format(map_result, "read1", 50);

  REQUIRE(results.size() == 1);
  // Path should be ">s1>s2>s3" (forward direction, original IDs)
  CHECK(results[0].graph_path == ">s1>s2>s3");
}

TEST_CASE("ResultFormatter handles primary and secondary mappings") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatter formatter(graph);

  piru::mapping::ReadMapResult map_result;

  // Primary mapping (higher score)
  piru::mapping::Mapping primary;
  primary.chain_score = 100.0;
  piru::mapping::ChainedAnchor a0;
  a0.node_id = 0; a0.offset = 0; a0.length = 4;
  a0.read_pos = 0;
  a0.path_id = 0;
  a0.ref_coord = 0;
  primary.anchors.push_back(a0);
  map_result.mappings.push_back(primary);

  // Secondary mapping (score >= 70% of primary to pass filter)
  piru::mapping::Mapping secondary;
  secondary.chain_score = 80.0;  // 80% of primary, passes 70% filter
  piru::mapping::ChainedAnchor a1;
  a1.node_id = 1; a1.offset = 0; a1.length = 4;
  a1.read_pos = 0;
  a1.path_id = 0;
  a1.ref_coord = 4;
  secondary.anchors.push_back(a1);
  map_result.mappings.push_back(secondary);

  auto results = formatter.format(map_result, "read1", 50);

  REQUIRE(results.size() == 2);

  // MAPQ equals raw chain score
  CHECK(results[0].mapq == 100);
  CHECK(results[1].mapq == 80);

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

TEST_CASE("ResultFormatter respects primary_only config") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatterConfig config;
  config.primary_only = true;
  piru::mapping::ResultFormatter formatter(graph, config);

  piru::mapping::ReadMapResult map_result;

  // Add two mappings
  piru::mapping::Mapping m1;
  m1.chain_score = 100.0;
  piru::mapping::ChainedAnchor a0;
  a0.node_id = 0; a0.offset = 0; a0.length = 4;
  a0.read_pos = 0;
  a0.path_id = 0;
  a0.ref_coord = 0;
  m1.anchors.push_back(a0);
  map_result.mappings.push_back(m1);

  piru::mapping::Mapping m2;
  m2.chain_score = 50.0;
  piru::mapping::ChainedAnchor a1;
  a1.node_id = 1; a1.offset = 0; a1.length = 4;
  a1.read_pos = 0;
  a1.path_id = 0;
  a1.ref_coord = 4;
  m2.anchors.push_back(a1);
  map_result.mappings.push_back(m2);

  auto results = formatter.format(map_result, "read1", 50);

  // Should only get primary
  CHECK(results.size() == 1);
}

TEST_CASE("ResultFormatter handles empty mappings") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatter formatter(graph);

  piru::mapping::ReadMapResult map_result;
  // No mappings

  auto results = formatter.format(map_result, "read1", 50);

  CHECK(results.empty());
}

TEST_CASE("ResultFormatter filters low-scoring secondaries") {
  auto graph = makeTestGraph();
  piru::mapping::ResultFormatterConfig config;
  config.min_secondary_ratio = 0.7;  // Filter secondaries below 70% of primary
  piru::mapping::ResultFormatter formatter(graph, config);

  piru::mapping::ReadMapResult map_result;

  // Primary mapping
  piru::mapping::Mapping primary;
  primary.chain_score = 100.0;
  piru::mapping::ChainedAnchor a0;
  a0.node_id = 0; a0.offset = 0; a0.length = 4;
  a0.read_pos = 0;
  a0.path_id = 0;
  a0.ref_coord = 0;
  primary.anchors.push_back(a0);
  map_result.mappings.push_back(primary);

  // Secondary that passes filter (70% of primary = 70)
  piru::mapping::Mapping secondary_good;
  secondary_good.chain_score = 75.0;  // Above threshold
  piru::mapping::ChainedAnchor a1;
  a1.node_id = 1; a1.offset = 0; a1.length = 4;
  a1.read_pos = 0;
  a1.path_id = 0;
  a1.ref_coord = 4;
  secondary_good.anchors.push_back(a1);
  map_result.mappings.push_back(secondary_good);

  // Secondary that fails filter
  piru::mapping::Mapping secondary_bad;
  secondary_bad.chain_score = 50.0;  // Below threshold
  piru::mapping::ChainedAnchor a2;
  a2.node_id = 2; a2.offset = 0; a2.length = 4;
  a2.read_pos = 0;
  a2.path_id = 0;
  a2.ref_coord = 8;
  secondary_bad.anchors.push_back(a2);
  map_result.mappings.push_back(secondary_bad);

  auto results = formatter.format(map_result, "read1", 50);

  // Should only get primary + good secondary (bad secondary filtered)
  CHECK(results.size() == 2);
  CHECK(results[0].mapq == 100);  // Primary chain score
  CHECK(results[1].mapq == 75);   // Good secondary chain score
}
