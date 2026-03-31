/**
 * test_bucket_indexer.cpp
 *
 * Microbenchmarks comparing bucket vs path-walk indexer seed output on
 * synthetic graphs covering different structural variant types.
 *
 * Each test builds a small FlatGraph with known structure, runs both
 * indexers, and compares the set of unique seed hashes produced.
 * Goal: understand which hashes differ and why.
 *
 * SPDX-License-Identifier: MIT
 */

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include "index/bucket_indexer.hpp"
#include "index/flat_graph.hpp"
#include "index/path_walk_indexer.hpp"
#include "index/simple_expand.hpp"
#include "io/graphs/graph.hpp"
#include "io/models/model.hpp"
#include "io/models/model_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"

using namespace piru;

namespace {

/* Helper: build an ImportedGraph from segments, edges, and paths. */
io::ImportedGraph makeGraph(
    const std::vector<std::pair<std::string, std::string>>& segments,
    const std::vector<std::tuple<std::string, std::string, bool, bool>>& edges,
    const std::vector<std::pair<std::string, std::vector<std::pair<std::string, bool>>>>& paths) {

    io::ImportedGraph g;
    for (const auto& [id, seq] : segments) {
        g.add_node({id, seq});
    }
    for (const auto& [from, to, from_rev, to_rev] : edges) {
        g.add_edge({from, to, from_rev, to_rev});
    }
    for (const auto& [name, steps] : paths) {
        io::ImportedPath p;
        p.name = name;
        for (const auto& [seg, rev] : steps) {
            p.steps.push_back({seg, rev});
        }
        g.add_path(std::move(p));
    }
    return g;
}

/* Helper: generate a deterministic DNA sequence of given length. */
std::string makeDNA(std::size_t len, int seed = 0) {
    static const char bases[] = "ACGT";
    std::string s(len, 'A');
    std::uint32_t state = static_cast<std::uint32_t>(seed) * 2654435761u + 1;
    for (std::size_t i = 0; i < len; ++i) {
        state = state * 1103515245u + 12345u;
        s[i] = bases[(state >> 16) & 3];
    }
    return s;
}

struct IndexResult {
    std::set<std::uint64_t> hashes;
    std::size_t total_seeds{0};
    std::size_t unique_hashes{0};
};

/* Run bucket indexer on an ImportedGraph and collect seed hashes. */
IndexResult runBucket(const io::ImportedGraph& imported, const io::KmerModel& model) {
    auto flat = index::simpleExpandFlat(imported);

    signal::FuzzyQuantizerConfig fq_cfg;
    fq_cfg.backend = "rh2";
    fq_cfg.pore_model = model.name();
    fq_cfg.diff = 0.35;  // match default RH2 diff filter
    auto fq = signal::make_fuzzy_quantizer(fq_cfg);

    signal::SeedExtractorConfig se_cfg;
    se_cfg.backend = "kmer";
    se_cfg.k = 8;
    se_cfg.stride = 1;
    se_cfg.qbits = 4;
    auto se = signal::make_seed_extractor(se_cfg);

    index::BucketIndexConfig cfg;
    cfg.seed_k = 8;
    cfg.seed_stride = 1;
    cfg.seed_freq_cutoff = 1.0;

    auto result = index::bucketIndex(flat, model, *fq, *se, cfg);

    IndexResult ir;
    // Collect hashes from FlatSeedStore via lookup iteration
    // We need to extract hashes -- use the raw accessors
    ir.unique_hashes = result.seed_store->size();
    ir.total_seeds = result.seeds_interior + result.seeds_boundary;

    // Collect all unique hashes by scanning the store's CSR
    const auto* flat_store = dynamic_cast<const index::FlatSeedStore*>(result.seed_store.get());
    if (flat_store) {
        for (const auto& h : flat_store->hashes()) {
            ir.hashes.insert(h);
        }
    }
    return ir;
}

/* Run path-walk indexer on an ImportedGraph and collect seed hashes. */
IndexResult runPathWalk(const io::ImportedGraph& imported, const io::KmerModel& model) {
    auto flat = index::simpleExpandFlat(imported);

    signal::FuzzyQuantizerConfig fq_cfg;
    fq_cfg.backend = "rh2";
    fq_cfg.pore_model = model.name();
    fq_cfg.diff = 0.0;
    auto fq = signal::make_fuzzy_quantizer(fq_cfg);

    signal::SeedExtractorConfig se_cfg;
    se_cfg.backend = "kmer";
    se_cfg.k = 8;
    se_cfg.stride = 1;
    se_cfg.qbits = 4;
    auto se = signal::make_seed_extractor(se_cfg);

    index::PathWalkIndexConfig cfg;
    cfg.seed_k = 8;
    cfg.seed_stride = 1;
    cfg.seed_freq_cutoff = 1.0;
    cfg.seed_freq_cap = 0.0;  // no subsampling

    auto result = index::pathWalkIndex(flat, model, *fq, *se, cfg);

    IndexResult ir;
    ir.unique_hashes = result.seeds_unique;
    ir.total_seeds = result.seeds_extracted;

    for (const auto& [hash, hits] : result.seed_store->data()) {
        ir.hashes.insert(hash);
    }
    return ir;
}

struct Comparison {
    std::size_t bucket_only{0};
    std::size_t pathwalk_only{0};
    std::size_t shared{0};
};

Comparison compareHashes(const IndexResult& bucket, const IndexResult& pathwalk) {
    Comparison c;
    for (auto h : bucket.hashes) {
        if (pathwalk.hashes.count(h)) c.shared++;
        else c.bucket_only++;
    }
    for (auto h : pathwalk.hashes) {
        if (!bucket.hashes.count(h)) c.pathwalk_only++;
    }
    return c;
}

void printComparison(const char* name, const IndexResult& b, const IndexResult& pw,
                     const Comparison& c) {
    MESSAGE(name << ":");
    MESSAGE("  bucket:    " << b.unique_hashes << " unique, " << b.total_seeds << " total");
    MESSAGE("  path-walk: " << pw.unique_hashes << " unique, " << pw.total_seeds << " total");
    MESSAGE("  shared: " << c.shared << "  bucket-only: " << c.bucket_only
            << "  pathwalk-only: " << c.pathwalk_only);
}

}  // namespace

// ============================================================================
// 1. LINEAR -- no variation, baseline
// ============================================================================
TEST_CASE("Bucket vs PathWalk: linear graph (no variation)") {
    auto model = io::load_builtin_model("r10.4");
    // Three nodes in a line, one path
    auto g = makeGraph(
        {{"s1", makeDNA(50, 1)}, {"s2", makeDNA(50, 2)}, {"s3", makeDNA(50, 3)}},
        {{"s1", "s2", false, false}, {"s2", "s3", false, false}},
        {{"path1", {{"s1", false}, {"s2", false}, {"s3", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("LINEAR", b, pw, c);

    // On a linear graph both should see the same hashes
    CHECK(c.pathwalk_only == 0);
    CHECK(c.bucket_only == 0);
}

// ============================================================================
// 2. SNP BUBBLE -- single base substitution
// ============================================================================
TEST_CASE("Bucket vs PathWalk: SNP bubble") {
    auto model = io::load_builtin_model("r10.4");
    // A -> B(1bp) -> D and A -> C(1bp) -> D
    auto g = makeGraph(
        {{"A", makeDNA(50, 10)}, {"B", "G"}, {"C", "T"}, {"D", makeDNA(50, 11)}},
        {{"A", "B", false, false}, {"A", "C", false, false},
         {"B", "D", false, false}, {"C", "D", false, false}},
        {{"p1", {{"A", false}, {"B", false}, {"D", false}}},
         {"p2", {{"A", false}, {"C", false}, {"D", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("SNP BUBBLE", b, pw, c);

    CHECK(c.shared > 0);
    // Log difference for investigation
    MESSAGE("  pathwalk-only hashes: " << c.pathwalk_only);
    MESSAGE("  bucket-only hashes: " << c.bucket_only);
}

// ============================================================================
// 3. SHORT INDEL -- 5bp insertion/deletion
// ============================================================================
TEST_CASE("Bucket vs PathWalk: short indel") {
    auto model = io::load_builtin_model("r10.4");
    // A -> B(5bp insertion) -> C, or A -> C directly
    auto g = makeGraph(
        {{"A", makeDNA(50, 20)}, {"B", "ACGTG"}, {"C", makeDNA(50, 21)}},
        {{"A", "B", false, false}, {"B", "C", false, false},
         {"A", "C", false, false}},
        {{"p_ins", {{"A", false}, {"B", false}, {"C", false}}},
         {"p_del", {{"A", false}, {"C", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("SHORT INDEL", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 4. STRUCTURAL VARIANT -- large insertion (500bp)
// ============================================================================
TEST_CASE("Bucket vs PathWalk: structural variant (large insertion)") {
    auto model = io::load_builtin_model("r10.4");
    auto g = makeGraph(
        {{"A", makeDNA(50, 30)}, {"B", makeDNA(500, 31)}, {"C", makeDNA(50, 32)}},
        {{"A", "B", false, false}, {"B", "C", false, false},
         {"A", "C", false, false}},
        {{"p_sv", {{"A", false}, {"B", false}, {"C", false}}},
         {"p_ref", {{"A", false}, {"C", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("SV LARGE INSERTION", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 5. CNV / TANDEM DUPLICATION -- path revisits a node
// ============================================================================
TEST_CASE("Bucket vs PathWalk: CNV (tandem duplication)") {
    auto model = io::load_builtin_model("r10.4");
    // Path visits B twice: A -> B -> B -> C
    auto g = makeGraph(
        {{"A", makeDNA(50, 40)}, {"B", makeDNA(50, 41)}, {"C", makeDNA(50, 42)}},
        {{"A", "B", false, false}, {"B", "B", false, false}, {"B", "C", false, false}},
        {{"p_cnv", {{"A", false}, {"B", false}, {"B", false}, {"C", false}}},
         {"p_ref", {{"A", false}, {"B", false}, {"C", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("CNV TANDEM DUP", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 6. INVERSION -- forward vs reverse of a segment
// ============================================================================
TEST_CASE("Bucket vs PathWalk: inversion") {
    auto model = io::load_builtin_model("r10.4");
    auto g = makeGraph(
        {{"A", makeDNA(50, 50)}, {"B", makeDNA(50, 51)}, {"C", makeDNA(50, 52)}},
        {{"A", "B", false, false}, {"A", "B", false, true},
         {"B", "C", false, false}, {"B", "C", true, false}},
        {{"p_fwd", {{"A", false}, {"B", false}, {"C", false}}},
         {"p_inv", {{"A", false}, {"B", true}, {"C", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("INVERSION", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 7. DISCONNECTED SUBGRAPHS -- two separate components
// ============================================================================
TEST_CASE("Bucket vs PathWalk: disconnected subgraphs") {
    auto model = io::load_builtin_model("r10.4");
    auto g = makeGraph(
        {{"A", makeDNA(50, 60)}, {"B", makeDNA(50, 61)},
         {"X", makeDNA(50, 62)}, {"Y", makeDNA(50, 63)}},
        {{"A", "B", false, false}, {"X", "Y", false, false}},
        {{"p1", {{"A", false}, {"B", false}}},
         {"p2", {{"X", false}, {"Y", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("DISCONNECTED", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 8. MULTI-BUBBLE -- chained SNPs
// ============================================================================
TEST_CASE("Bucket vs PathWalk: multi-bubble (chained SNPs)") {
    auto model = io::load_builtin_model("r10.4");
    // A -> B/C -> D -> E/F -> G
    auto g = makeGraph(
        {{"A", makeDNA(50, 70)}, {"B", "G"}, {"C", "T"},
         {"D", makeDNA(50, 71)}, {"E", "A"}, {"F", "C"},
         {"G", makeDNA(50, 72)}},
        {{"A", "B", false, false}, {"A", "C", false, false},
         {"B", "D", false, false}, {"C", "D", false, false},
         {"D", "E", false, false}, {"D", "F", false, false},
         {"E", "G", false, false}, {"F", "G", false, false}},
        {{"p1", {{"A", false}, {"B", false}, {"D", false}, {"E", false}, {"G", false}}},
         {"p2", {{"A", false}, {"C", false}, {"D", false}, {"F", false}, {"G", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("MULTI-BUBBLE", b, pw, c);
    MESSAGE("  pathwalk-only: " << c.pathwalk_only << "  bucket-only: " << c.bucket_only);
}

// ============================================================================
// 9. TINY NODE -- shorter than pore_k, all seeds from boundary
// ============================================================================
TEST_CASE("Bucket vs PathWalk: tiny node (shorter than pore_k)") {
    auto model = io::load_builtin_model("r10.4");
    // B is only 3bp -- too short for any interior seeds with k=9
    auto g = makeGraph(
        {{"A", makeDNA(50, 80)}, {"B", "ACG"}, {"C", makeDNA(50, 81)}},
        {{"A", "B", false, false}, {"B", "C", false, false}},
        {{"p1", {{"A", false}, {"B", false}, {"C", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("TINY NODE", b, pw, c);
    CHECK(c.pathwalk_only == 0);
    CHECK(c.bucket_only == 0);
}

// ============================================================================
// 10. LONG SINGLE NODE -- no boundaries needed
// ============================================================================
TEST_CASE("Bucket vs PathWalk: long single node") {
    auto model = io::load_builtin_model("r10.4");
    auto g = makeGraph(
        {{"mega", makeDNA(1000, 90)}},
        {},
        {{"p1", {{"mega", false}}}});

    auto b = runBucket(g, *model);
    auto pw = runPathWalk(g, *model);
    auto c = compareHashes(b, pw);
    printComparison("LONG SINGLE NODE", b, pw, c);

    // Single node, no boundaries -- must be identical
    CHECK(c.pathwalk_only == 0);
    CHECK(c.bucket_only == 0);
}
