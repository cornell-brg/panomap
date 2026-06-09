// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

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
