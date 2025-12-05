// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"

#include <doctest/doctest.h>

using namespace piru::signal;

TEST_CASE("K-mer seed extractor emits sliding window hashes") {
    FuzzyQuantizedSignal quantized;
    quantized.tokens = {1, 2, 3};

    auto extractor = make_seed_extractor(
        SeedExtractorConfig{.backend = "kmer", .k = 2, .stride = 1, .qbits = 4});
    auto seeds = extractor->extract(quantized, nullptr);

    REQUIRE(seeds.seeds.size() == 2);
    CHECK(seeds.seeds[0].hash == 3911557750ULL);  // hash([1,2])
    CHECK(seeds.seeds[0].position == 0);
    CHECK(seeds.seeds[0].span == 2);
    CHECK(seeds.seeds[1].hash == 447855090ULL);  // hash([2,3])
    CHECK(seeds.seeds[1].position == 1);
    CHECK(seeds.seeds[1].span == 2);
}

