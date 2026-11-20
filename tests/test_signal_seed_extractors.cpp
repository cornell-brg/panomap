// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "signal/seed_extractors/seed_hash_util.hpp"

#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"

using namespace panomap::signal;

TEST_CASE("K-mer seed extractor emits sliding window hashes") {
  TokenizedSignal quantized;
  quantized.tokens = {1, 2, 3};

  auto extractor =
      make_seed_extractor(SeedExtractorConfig{.backend = "kmer", .k = 2, .stride = 1, .qbits = 4});
  auto seeds = extractor->extract(quantized);

  REQUIRE(seeds.seeds.size() == 2);
  CHECK(seeds.seeds[0].hash == 3911557750ULL);  // hash([1,2])
  CHECK(seeds.seeds[0].position == 0);
  CHECK(seeds.seeds[0].length == 2);           // Initialized to k
  CHECK(seeds.seeds[1].hash == 447855090ULL);  // hash([2,3])
  CHECK(seeds.seeds[1].position == 1);
  CHECK(seeds.seeds[1].length == 2);  // Initialized to k
}

TEST_CASE("Seed struct length field is initialized correctly") {
  TokenizedSignal quantized;
  quantized.tokens = {1, 2, 3, 4, 5, 6};

  auto extractor =
      make_seed_extractor(SeedExtractorConfig{.backend = "kmer", .k = 3, .stride = 2, .qbits = 4});
  auto seeds = extractor->extract(quantized);

  REQUIRE(seeds.seeds.size() == 2);

  // All seeds should have length initialized to k
  for (const auto& seed : seeds.seeds) {
    CHECK(seed.length == 3);
  }
}

TEST_CASE("Seed length can be updated after merging") {
  // Test that length can be modified (e.g., after merging multiple seeds)
  Seed seed;
  seed.hash = 12345;
  seed.position = 10;
  seed.length = 20;  // Initially set to k

  CHECK(seed.length == 20);

  // Simulate merging: length increases to cover merged region
  seed.length = 35;
  CHECK(seed.length == 35);  // Updated coverage after hypothetical merge
}

TEST_CASE("Minimizer seed extractor emits every position tied for the window minimum") {
  TokenizedSignal quantized;
  quantized.tokens = {5, 5, 5, 5};  // three equal 2-mers in one window of 3

  auto extractor = make_seed_extractor(
      SeedExtractorConfig{.backend = "minimizer", .k = 2, .qbits = 4, .window = 3});
  auto seeds = extractor->extract(quantized);

  REQUIRE(seeds.seeds.size() == 3);
  CHECK(seeds.seeds[0].position == 0);
  CHECK(seeds.seeds[1].position == 1);
  CHECK(seeds.seeds[2].position == 2);
  CHECK(seeds.seeds[0].hash == seeds.seeds[2].hash);
}

TEST_CASE("seed hash: every token of a seed longer than 8 tokens counts") {
  using panomap::signal::hashKmer;
  using panomap::signal::makeMask;
  const std::uint32_t qbits = 4;
  const std::uint64_t token_mask = makeMask(qbits);
  // Two 10-token seeds that differ only in their first token.
  const std::int16_t a[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const std::int16_t b[10] = {11, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const std::uint64_t mask10 = makeMask(10 * qbits);
  CHECK(hashKmer(a, 0, 10, qbits, token_mask, mask10, true) !=
        hashKmer(b, 0, 10, qbits, token_mask, mask10, true));
  // An 8-token seed keeps its 32-bit hash.
  const std::uint64_t mask8 = makeMask(8 * qbits);
  CHECK(hashKmer(a, 0, 8, qbits, token_mask, mask8, true) <= 0xFFFFFFFFULL);
}
