// SPDX-License-Identifier: MIT

#include "index/seed_builder.hpp"
#include "index/transform_dbg.hpp"
#include "index/squigglize.hpp"
#include "index/pseudo_linearize.hpp"

#include <doctest/doctest.h>

#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"

using namespace piru;

namespace {

class ConstKmerModel : public io::KmerModel {
public:
    explicit ConstKmerModel(int k, double mean) : k_(k), mean_(mean) {}

    std::string name() const override { return "const"; }
    int k() const override { return k_; }
    bool lookup(const std::string&, double& mean) const override {
        mean = mean_;
        return true;
    }

private:
    int k_;
    double mean_;
};

index::AlnGraph makeGraph(const std::vector<std::string>& seqs) {
    index::AlnGraph g;
    for (const auto& seq : seqs) {
        index::AlnNode n;
        n.label = seq;
        n.original_id = seq;
        n.sequence = seq;
        g.addNode(n);
    }
    return g;
}

}  // namespace

TEST_CASE("SCC detection marks cycles and orders components") {
    index::AlnGraph g;
    g.addNode({}); g.addNode({}); g.addNode({}); g.addNode({}); g.addNode({});
    g.addEdge({0, 1, 0});
    g.addEdge({1, 2, 0});
    g.addEdge({2, 0, 0});  // cycle among 0,1,2
    g.addEdge({3, 0, 0});  // enters cycle
    g.addEdge({2, 4, 0});  // exits cycle

    const auto scc = index::computeScc(g);
    REQUIRE(scc.component.size() == 5);
    REQUIRE(scc.in_scc.size() == 5);

    // Cycle nodes marked and share component.
    CHECK(scc.in_scc[0] == true);
    CHECK(scc.in_scc[1] == true);
    CHECK(scc.in_scc[2] == true);
    CHECK(scc.component[0] == scc.component[1]);
    CHECK(scc.component[1] == scc.component[2]);

    // Node 3 upstream of cycle, node 4 downstream.
    CHECK(scc.component[3] < scc.component[0]);
    CHECK(scc.component[2] < scc.component[4]);
}

TEST_CASE("Tip folding marks tips and keeps main chain intact") {
    // Components: {4,5} tip, {0,1,2} main, {6} downstream target of main.
    index::AlnGraph g;
    g.addNode({}); g.addNode({}); g.addNode({}); g.addNode({});
    g.addNode({}); g.addNode({}); g.addNode({});
    g.addEdge({0, 1, 0});
    g.addEdge({1, 2, 0});
    g.addEdge({2, 6, 0});  // main chain has outgoing edge to non-tip component

    index::SccResult scc;
    scc.component = {1, 1, 1, 2, 0, 0, 2};
    scc.in_scc = {false, false, false, false, false, false, false};
    scc.components = 3;

    auto tip = index::chainTips(g, scc);
    CHECK(tip.ignorable_tip[4] == true);
    CHECK(tip.ignorable_tip[5] == true);
    CHECK(tip.ignorable_tip[0] == false);
    CHECK(tip.ignorable_tip[1] == false);
    CHECK(tip.ignorable_tip[2] == false);
}

TEST_CASE("Cycle folding collapses simple 2-cycles") {
    index::AlnGraph g;
    g.addNode({}); g.addNode({}); g.addEdge({0, 1, 0}); g.addEdge({1, 0, 0});

    index::TipFoldingResult tip(2);
    tip.ignorable_tip.assign(2, false);

    index::chainCycles(g, tip);
    CHECK(tip.uf.find(0) == tip.uf.find(1));
    CHECK((tip.ignorable_tip[0] || tip.ignorable_tip[1]) == true);
}

TEST_CASE("Superbubble chaining merges start/end and interior") {
    // Graph: 0 -> 1 -> 3
    //        0 -> 2 -> 3
    // Superbubble between 0 (start) and 3 (exit).
    index::AlnGraph g;
    g.addNode({}); g.addNode({}); g.addNode({}); g.addNode({});
    g.addEdge({0, 1, 0});
    g.addEdge({0, 2, 0});
    g.addEdge({1, 3, 0});
    g.addEdge({2, 3, 0});

    index::SccResult scc;
    scc.component = {0, 1, 1, 2};
    scc.in_scc = {false, false, false, false};
    scc.components = 3;

    index::TipFoldingResult tips(4);
    tips.ignorable_tip.assign(4, false);

    auto sb = index::chainSuperbubbles(g, scc, tips);
    const auto rep0 = sb.uf.find(0);
    CHECK(sb.uf.find(3) == rep0);  // start and end merged
    CHECK(sb.uf.find(1) == rep0);  // interior merged
    CHECK(sb.uf.find(2) == rep0);  // interior merged
}

TEST_CASE("DBG transform splits nodes and trims overlaps") {
    // Simple DBG with k_graph=5, pore_k=3 => k_delta=2 => trim last 2 bases
    io::ImportedGraph imported;
    imported.add_node({.id = "n1", .sequence = "ACGTA"});
    imported.add_node({.id = "n2", .sequence = "CGTAC"});
    imported.add_edge({.from = "n1",
                       .to = "n2",
                       .from_reverse = false,
                       .to_reverse = false,
                       .overlap = "4",
                       .overlap_bases = 4});

    const std::size_t graph_k = 5;
    const std::size_t pore_k = 3;  // keep 2bp overlap

    auto aln = index::transformDbg(imported, graph_k, pore_k);

    // Two nodes become four directional nodes.
    REQUIRE(aln.nodeCount() == 4);

    // Sequences trimmed to length 3 (ACGTA -> ACG, CGTAC -> CGT).
    CHECK(aln.node(0).sequence == "ACG");      // n1 forward
    CHECK(aln.node(1).sequence == "CGT");      // n1 reverse (revcomp of ACG)
    CHECK(aln.node(2).sequence == "CGT");      // n2 forward
    CHECK(aln.node(3).sequence == "ACG");      // n2 reverse (revcomp of CGT)

    // One edge mapped forward with overlap reduced from 4 -> 2.
    CHECK(aln.edgeCount() == 1);
    CHECK(aln.outgoing(0).size() == 1);
    const auto succ = aln.outgoing(0).front();
    CHECK(succ == 2);  // n1 fwd -> n2 fwd
}

TEST_CASE("squigglize produces per-node signals with expected lengths") {
    // Model: k=1, constant mean 5.0 → raw signals all 5.0 → normalized to zero.
    ConstKmerModel model(1, 5.0);
    signal::FuzzyQuantizerConfig fq_cfg;
    fq_cfg.backend = "passthrough";
    auto fuzzy = signal::make_fuzzy_quantizer(fq_cfg);

    signal::AlignmentQuantizerConfig aq_cfg;
    aq_cfg.backend = "passthrough";
    auto align_q = signal::make_alignment_quantizer(aq_cfg);

    auto graph = makeGraph({"AAAAA", "AAA"});

    const auto result =
        index::squigglizeAndQuantize(graph, model, *fuzzy, *align_q);

    REQUIRE(result.fuzzy_signals.size() == 2);
    REQUIRE(result.alignment_signals.size() == 2);

    // k=1 → length = seq_len - k + 1.
    CHECK(result.fuzzy_signals[0].tokens.size() == 5);
    CHECK(result.fuzzy_signals[1].tokens.size() == 3);

    // Constant model + global normalization → zeros after quantization.
    for (auto t : result.fuzzy_signals[0].tokens) CHECK(t == 0);
    for (auto t : result.fuzzy_signals[1].tokens) CHECK(t == 0);
}

TEST_CASE("seed builder computes frequencies and threshold") {
    // Build two fuzzy signals with simple tokens.
    signal::FuzzyQuantizedSignal s0;
    s0.tokens = {1, 1, 1, 1};  // one hash frequency=4 (with k=1, stride=1)
    signal::FuzzyQuantizedSignal s1;
    s1.tokens = {2, 3};        // two hashes frequency=1 each

    std::vector<signal::FuzzyQuantizedSignal> signals{s0, s1};

    signal::SeedExtractorConfig cfg;
    cfg.backend = "kmer";
    cfg.k = 1;
    cfg.stride = 1;
    cfg.qbits = 16;
    auto extractor = signal::make_seed_extractor(cfg);

    index::SeedBuildConfig build_cfg;
    build_cfg.keep_least_frequent_fraction = 0.5;  // keep up to 50th percentile

    const auto store = index::buildSeedStore(signals, *extractor, build_cfg);

    CHECK(store.size() == 3);                // hashes for 1, 2, 3
    CHECK(store.max_hash_frequency() == 4);  // hash(1) seen 4 times
    CHECK(store.frequency_threshold() == 2); // 50th percentile (1,1,4) → 1 + 1
}
