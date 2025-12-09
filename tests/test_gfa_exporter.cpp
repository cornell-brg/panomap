#include "doctest/doctest.h"

#include <fstream>
#include <string>
#include <vector>

#include "io/graphs/gfa_exporter.hpp"
#include "index/aln_graph.hpp"
#include "io/graphs/graph.hpp"

using namespace piru;

TEST_CASE("GfaExporter dumpImportedGraph") {
    io::ImportedGraph graph;
    graph.add_node({"n1", "ACGT"});
    graph.add_node({"n2", "TGCA"});
    graph.add_edge({"n1", "n2", false, false, "4M", 4});
    io::ImportedPath path;
    path.name = "p1";
    path.steps.push_back({"n1", false});
    path.steps.push_back({"n2", false});
    path.overlaps.push_back("4M");
    graph.add_path(path);

    const std::string test_file = "test_imported.gfa";
    GfaExporter::dumpImportedGraph(graph, test_file);

    std::ifstream in(test_file);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(in, line)) {
        lines.push_back(line);
    }

    CHECK(lines.size() == 5);
    CHECK(lines[0] == "H\tVN:Z:1.0");
    CHECK(lines[1] == "S\tn1\tACGT\tLN:i:4");
    CHECK(lines[2] == "S\tn2\tTGCA\tLN:i:4");
    CHECK(lines[3] == "L\tn1\t+\tn2\t+\t4M");
    CHECK(lines[4] == "P\tp1\tn1+,n2+\t4M");
}

TEST_CASE("GfaExporter dumpAlnGraph") {
    index::AlnGraph graph;
    index::AlnNode node1;
    node1.id = 0;
    node1.label = "n1_F";
    node1.sequence = "ACGT";
    node1.chain_id = 1;
    node1.linear_position = 100;
    graph.addNode(node1);

    index::AlnNode node2;
    node2.id = 1;
    node2.label = "n2_F";
    node2.sequence = "TGCA";
    node2.chain_id = 1;
    node2.linear_position = 200;
    graph.addNode(node2);
    
    graph.addEdge({0, 1, 3});

    index::AlnPath path;
    path.name = "path1";
    path.steps.push_back({"n1_F", false});
    path.steps.push_back({"n2_F", false});
    path.overlaps.push_back(3);
    graph.addPath(path);

    const std::string test_file = "test_aln.gfa";
    GfaExporter::dumpAlnGraph(graph, test_file);

    std::ifstream in(test_file);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(in, line)) {
        lines.push_back(line);
    }
    
    CHECK(lines.size() == 5);
    CHECK(lines[0] == "H\tVN:Z:1.0");
    CHECK(lines[1] == "S\tn1_F\tACGT\tLN:i:4\tci:i:1\tlc:i:100");
    CHECK(lines[2] == "S\tn2_F\tTGCA\tLN:i:4\tci:i:1\tlc:i:200");
    CHECK(lines[3] == "L\tn1_F\t+\tn2_F\t+\t3M");
    CHECK(lines[4] == "P\tpath1\tn1_F+,n2_F+\t3M");
}

TEST_CASE("GfaExporter dumpAlnGraph RawSignal") {
    index::AlnGraph graph;
    index::AlnNode node1;
    node1.id = 0;
    node1.label = "n1_F";
    graph.addNode(node1);

    std::vector<std::vector<float>> signals = {{1.1, 2.2, 3.3}};
    const std::string test_file = "test_aln_raw.gfa";
    GfaExporter::dumpAlnGraph(graph, test_file, AlnGraphDumpMode::RawSignal, &signals);

    std::ifstream in(test_file);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(in, line)) {
        lines.push_back(line);
    }

    CHECK(lines.size() == 2);
    CHECK(lines[0] == "H\tVN:Z:1.0");
    CHECK(lines[1] == "S\tn1_F\t1.1,2.2,3.3\tLN:i:3\tst:Z:raw_signal");
}

TEST_CASE("GfaExporter dumpAlnGraph FuzzyQuantized") {
    index::AlnGraph graph;
    index::AlnNode node1;
    node1.id = 0;
    node1.label = "n1_F";
    graph.addNode(node1);

    std::vector<std::vector<int16_t>> signals = {{10, 20, -30}};
    const std::string test_file = "test_aln_fuzzy.gfa";
    GfaExporter::dumpAlnGraph(graph, test_file, AlnGraphDumpMode::FuzzyQuantized, &signals);

    std::ifstream in(test_file);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(in, line)) {
        lines.push_back(line);
    }

    CHECK(lines.size() == 2);
    CHECK(lines[0] == "H\tVN:Z:1.0");
    CHECK(lines[1] == "S\tn1_F\t10,20,-30\tLN:i:3\tst:Z:fuzzy_quant");
}

TEST_CASE("GfaExporter dumpAlnGraph AlnQuantized") {
    index::AlnGraph graph;
    index::AlnNode node1;
    node1.id = 0;
    node1.label = "n1_F";
    graph.addNode(node1);

    std::vector<std::vector<int16_t>> signals = {{100, 200, -300}};
    const std::string test_file = "test_aln_aln.gfa";
    GfaExporter::dumpAlnGraph(graph, test_file, AlnGraphDumpMode::AlnQuantized, &signals);

    std::ifstream in(test_file);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(in, line)) {
        lines.push_back(line);
    }

    CHECK(lines.size() == 2);
    CHECK(lines[0] == "H\tVN:Z:1.0");
    CHECK(lines[1] == "S\tn1_F\t100,200,-300\tLN:i:3\tst:Z:aln_quant");
}
