// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>
#include <sstream>
#include <string>

#include "io/reads/read_provider_factory.hpp"
#include "mapping/batch_mapper.hpp"

using namespace piru;

TEST_CASE("Batch mapper processes blow5 reads and emits seeds") {
#ifdef PIRU_HAS_SLOW5
  const std::string path = "tests/data/HLA/test_reads/quick_r9_1k.blow5";
  auto provider = io::make_read_provider(path);
  REQUIRE(provider != nullptr);

  mapping::BatchMapperConfig config;
  config.num_threads = 1;

  std::ostringstream sink;
  mapping::BatchMapper mapper(*provider, config, sink);
  const auto stats = mapper.process_all();

  CHECK(stats.reads_processed == 5);
  CHECK(stats.batches == 1);

  const auto output = sink.str();
  CHECK_FALSE(output.empty());
  CHECK(output.find("seeds=") != std::string::npos);
#else
  MESSAGE("PIRU_HAS_SLOW5 not set; skipping batch mapper integration test");
  CHECK(true);
#endif
}
