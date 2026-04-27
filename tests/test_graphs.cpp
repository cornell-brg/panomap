#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

#include "core/io/graphs/graph.hpp"
#include "core/io/graphs/graph_loader_factory.hpp"

TEST_CASE("GFA loader parses segments and links") {
  const auto src_dir = std::filesystem::path(__FILE__).parent_path();
  const auto gfa_path = src_dir / "data/graphs/sample.gfa";

  auto loader = piru::io::make_graph_loader(gfa_path.string());
  REQUIRE(loader != nullptr);
  CHECK(loader->get_format_name() == "gfa");

  piru::io::ImportedGraph graph;
  REQUIRE(loader->load(graph));
  CHECK(graph.nodes.size() == 2);
  CHECK(graph.edges.size() == 1);
  CHECK(graph.paths.size() == 1);

  const auto* first = graph.find_node("1");
  REQUIRE(first != nullptr);
  CHECK(first->sequence == "ACGT");

  const auto* second = graph.find_node("2");
  REQUIRE(second != nullptr);
  CHECK(second->sequence == "GG");

  const auto& edge = graph.edges.front();
  CHECK(edge.from == "1");
  CHECK(edge.to == "2");
  CHECK_FALSE(edge.from_reverse);
  CHECK(edge.to_reverse);
  CHECK(edge.overlap_bases.has_value());
  CHECK(*edge.overlap_bases == 0);

  const auto& imported_path = graph.paths.front();
  CHECK(imported_path.name == "path1");
  REQUIRE(imported_path.steps.size() == 2);
  CHECK(imported_path.steps[0].segment_id == "1");
  CHECK_FALSE(imported_path.steps[0].is_reverse);
  CHECK(imported_path.steps[1].segment_id == "2");
  CHECK(imported_path.steps[1].is_reverse);
  REQUIRE(imported_path.overlaps.size() == 1);
  // CHECK(imported_path.overlaps.front().has_value());
  // CHECK(*imported_path.overlaps.front() == 0);
}

TEST_CASE("ImportedGraph basic helpers") {
  piru::io::ImportedGraph graph;

  graph.add_node({"id1", "AC"});
  graph.add_edge({"id1", "id2", false, true, "5M", 5});
  graph.add_path({"p", {{"id1", false}}, {}});

  const auto* found = graph.find_node("id1");
  REQUIRE(found != nullptr);
  CHECK(found->sequence == "AC");
  CHECK(graph.edges.front().overlap == "5M");
  CHECK(graph.edges.front().overlap_bases.has_value());
  CHECK(*graph.edges.front().overlap_bases == 5);
}

