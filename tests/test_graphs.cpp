#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>

#include "io/graphs/graph.hpp"
#include "io/graphs/graph_loader_factory.hpp"

#ifdef PIRU_HAS_LIBVGIO
#include <vg/io/stream.hpp>
#include <vg/vg.pb.h>
#endif

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
    CHECK(graph.flavor == piru::io::ImportedGraphFlavor::kUnknown);

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

TEST_CASE("vg loader parses binary vg via libvgio") {
#ifdef PIRU_HAS_LIBVGIO
    vg::Graph proto_graph;
    auto* node1 = proto_graph.add_node();
    node1->set_id(1);
    node1->set_sequence("AC");
    auto* node2 = proto_graph.add_node();
    node2->set_id(2);
    node2->set_sequence("T");
    auto* edge = proto_graph.add_edge();
    edge->set_from(1);
    edge->set_to(2);
    edge->set_from_start(false);
    edge->set_to_end(true);
    auto* path = proto_graph.add_path();
    path->set_name("p1");
    auto* mapping1 = path->add_mapping();
    mapping1->set_rank(1);
    auto* pos1 = mapping1->mutable_position();
    pos1->set_node_id(1);
    pos1->set_is_reverse(false);
    auto* mapping2 = path->add_mapping();
    mapping2->set_rank(2);
    auto* pos2 = mapping2->mutable_position();
    pos2->set_node_id(2);
    pos2->set_is_reverse(false);

    const auto tmp_path =
        std::filesystem::temp_directory_path() / "piru_test_sample.vg";
    std::ofstream out(tmp_path, std::ios::binary);
    vg::io::write<vg::Graph>(out, 1, [&proto_graph](std::size_t) { return proto_graph; }, false);
    vg::io::finish(out, false);
    out.close();

    auto loader = piru::io::make_graph_loader(tmp_path.string());
    REQUIRE(loader != nullptr);
    CHECK(loader->get_format_name() == "vg");

    piru::io::ImportedGraph graph;
    REQUIRE(loader->load(graph));
    CHECK(graph.nodes.size() == 2);
    CHECK(graph.edges.size() == 1);
    CHECK(graph.paths.size() == 1);

    const auto* first = graph.find_node("1");
    REQUIRE(first != nullptr);
    CHECK(first->sequence == "AC");

    const auto* second = graph.find_node("2");
    REQUIRE(second != nullptr);
    CHECK(second->sequence == "T");

    const auto& parsed_edge = graph.edges.front();
    CHECK(parsed_edge.from == "1");
    CHECK(parsed_edge.to == "2");
    CHECK_FALSE(parsed_edge.from_reverse);
    CHECK(parsed_edge.to_reverse);
    CHECK(parsed_edge.overlap == "0");
    CHECK_FALSE(parsed_edge.overlap_bases.has_value());
    const auto& parsed_path = graph.paths.front();
    CHECK(parsed_path.name == "p1");
    CHECK(parsed_path.steps.size() == 2);

#else
    MESSAGE("PIRU_HAS_LIBVGIO not set; skipping vg loader test");
    CHECK(true);
#endif
}
