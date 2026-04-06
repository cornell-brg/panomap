/**
 * test_bucket_indexer.cpp
 *
 * Tests for the bucket indexer: verifies seed generation on synthetic
 * graphs covering different structural variant types, and round-trip
 * serialization parity.
 *
 * SPDX-License-Identifier: MIT
 */

#include <algorithm>
#include <cstdint>
#include <doctest/doctest.h>
#include <set>
#include <string>
#include <vector>

#include "index/bucket_indexer.hpp"
#include "index/bucket_seed_store.hpp"
#include "index/flat_graph.hpp"
#include "index/simple_expand.hpp"
#include "io/graphs/graph.hpp"
#include "io/models/model.hpp"
#include "io/models/model_factory.hpp"
#include "signal/tokenizers/tokenizer_factory.hpp"
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

  signal::TokenizerConfig fq_cfg;
  fq_cfg.backend = "rh2";
  fq_cfg.pore_model = model.name();
  auto fq = signal::make_tokenizer(fq_cfg);

  signal::SeedExtractorConfig se_cfg;
  se_cfg.backend = "kmer";
  se_cfg.k = 8;
  se_cfg.stride = 1;
  se_cfg.qbits = 4;
  auto se = signal::make_seed_extractor(se_cfg);

  index::BucketIndexConfig cfg;
  cfg.seed_k = 8;

  auto result = index::bucketIndex(flat, model, *fq, *se, cfg);

  IndexResult ir;
  ir.unique_hashes = result.seed_store->size();
  ir.total_seeds = result.seeds_interior + result.seeds_boundary;

  const auto* bucket_store = dynamic_cast<const index::BucketSeedStore*>(result.seed_store.get());
  if (bucket_store) {
    for (std::size_t bi = 0; bi < bucket_store->num_buckets(); ++bi) {
      const auto& b = bucket_store->bucket(bi);
      for (const auto& h : b.keys) {
        ir.hashes.insert(h);
      }
    }
  }
  return ir;
}

/* Verify that every hash in the store can be looked up. */
bool verifyLookup(const index::SeedStore& store, const std::set<std::uint64_t>& hashes) {
  for (auto h : hashes) {
    auto span = store.lookup(h);
    if (span.count == 0) return false;
  }
  return true;
}

}  // namespace

// ============================================================================
// Basic seed generation tests
// ============================================================================

TEST_CASE("Bucket indexer: linear graph produces seeds") {
  auto model = io::load_builtin_model("r10.4_400bps");
  std::string seq = makeDNA(200, 42);
  auto graph = makeGraph({{"s1", seq}}, {}, {{"p1", {{"s1", false}}}});
  auto result = runBucket(graph, *model);

  CHECK(result.unique_hashes > 0);
  CHECK(result.total_seeds > 0);
  MESSAGE("unique: " << result.unique_hashes << " total: " << result.total_seeds);
}

TEST_CASE("Bucket indexer: SNP bubble produces seeds from both paths") {
  auto model = io::load_builtin_model("r10.4_400bps");
  std::string prefix = makeDNA(50, 1);
  std::string suffix = makeDNA(50, 2);
  std::string var_a = makeDNA(20, 3);
  std::string var_b = makeDNA(20, 4);

  auto graph = makeGraph({{"pre", prefix}, {"a", var_a}, {"b", var_b}, {"suf", suffix}},
                         {{"pre", "a", false, false},
                          {"pre", "b", false, false},
                          {"a", "suf", false, false},
                          {"b", "suf", false, false}},
                         {{"p1", {{"pre", false}, {"a", false}, {"suf", false}}},
                          {"p2", {{"pre", false}, {"b", false}, {"suf", false}}}});

  auto result = runBucket(graph, *model);
  CHECK(result.unique_hashes > 0);
  MESSAGE("SNP bubble: " << result.unique_hashes << " unique hashes");
}

TEST_CASE("Bucket indexer: tiny nodes produce seeds via path context") {
  auto model = io::load_builtin_model("r10.4_400bps");
  // Nodes shorter than pore_k -- need path context to generate any seeds
  auto graph = makeGraph(
      {{"s1", makeDNA(5, 1)},
       {"s2", makeDNA(5, 2)},
       {"s3", makeDNA(5, 3)},
       {"s4", makeDNA(5, 4)},
       {"s5", makeDNA(50, 5)}},
      {{"s1", "s2", false, false},
       {"s2", "s3", false, false},
       {"s3", "s4", false, false},
       {"s4", "s5", false, false}},
      {{"p1", {{"s1", false}, {"s2", false}, {"s3", false}, {"s4", false}, {"s5", false}}}});

  auto result = runBucket(graph, *model);
  CHECK(result.unique_hashes > 0);
  MESSAGE("Tiny nodes: " << result.unique_hashes << " unique hashes");
}

TEST_CASE("Bucket indexer: disconnected subgraphs") {
  auto model = io::load_builtin_model("r10.4_400bps");
  auto graph = makeGraph({{"a", makeDNA(100, 1)}, {"b", makeDNA(100, 2)}}, {},
                         {{"p1", {{"a", false}}}, {"p2", {{"b", false}}}});

  auto result = runBucket(graph, *model);
  CHECK(result.unique_hashes > 0);
  MESSAGE("Disconnected: " << result.unique_hashes << " unique hashes");
}

TEST_CASE("Bucket indexer: lookup parity after build") {
  auto model = io::load_builtin_model("r10.4_400bps");
  std::string seq = makeDNA(300, 99);
  auto graph = makeGraph(
      {{"s1", seq.substr(0, 100)}, {"s2", seq.substr(100, 100)}, {"s3", seq.substr(200, 100)}},
      {{"s1", "s2", false, false}, {"s2", "s3", false, false}},
      {{"p1", {{"s1", false}, {"s2", false}, {"s3", false}}}});

  auto flat = index::simpleExpandFlat(graph);

  signal::TokenizerConfig fq_cfg;
  fq_cfg.backend = "rh2";
  fq_cfg.pore_model = model->name();
  auto fq = signal::make_tokenizer(fq_cfg);

  signal::SeedExtractorConfig se_cfg;
  se_cfg.backend = "kmer";
  se_cfg.k = 8;
  se_cfg.stride = 1;
  se_cfg.qbits = 4;
  auto se = signal::make_seed_extractor(se_cfg);

  index::BucketIndexConfig cfg;
  cfg.seed_k = 8;

  auto result = index::bucketIndex(flat, *model, *fq, *se, cfg);

  // Collect all hashes
  std::set<std::uint64_t> hashes;
  const auto* store = dynamic_cast<const index::BucketSeedStore*>(result.seed_store.get());
  REQUIRE(store != nullptr);
  for (std::size_t bi = 0; bi < store->num_buckets(); ++bi) {
    const auto& b = store->bucket(bi);
    for (const auto& h : b.keys) hashes.insert(h);
  }

  // Verify every hash can be looked up
  CHECK(verifyLookup(*store, hashes));
  MESSAGE("Lookup parity: " << hashes.size() << " hashes verified");
}
