// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include <cstdint>
#include <limits>
#include <vector>

#include "signal/tokenizers/tokenizer_factory.hpp"

using namespace panomap::signal;

TEST_CASE("RH2 tokenizer bins within fine range") {
  NormalizedSignal norm;
  norm.samples = {0.0f};

  TokenizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -2.0f;
  cfg.fine_max = 2.0f;
  cfg.fine_range = 0.4f;
  cfg.qbits = 4;

  auto tokenizer = make_tokenizer(cfg);
  auto out = tokenizer->quantize(norm);

  REQUIRE(out.tokens.size() == 1);
  CHECK(out.tokens[0] == 3);  // rawhash2 formula -> bucket ~3 for 0.0
}

TEST_CASE("RH2 tokenizer clamps and maps extremes") {
  NormalizedSignal norm;
  norm.samples = {-3.0f, 3.0f};

  TokenizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -2.0f;
  cfg.fine_max = 2.0f;
  cfg.fine_range = 0.4f;
  cfg.qbits = 4;

  auto tokenizer = make_tokenizer(cfg);
  auto out = tokenizer->quantize(norm);

  REQUIRE(out.tokens.size() == 2);
  CHECK(out.tokens[0] == 6);   // lower extreme maps to coarse bucket
  CHECK(out.tokens[1] == 15);  // upper extreme maps to top bucket
}

TEST_CASE("RH2 tokenizer collapses repeated tokens when asked") {
  NormalizedSignal norm;
  // 0.0 and 0.01 fall in the same bucket; 3.0 in another; NaN is a sentinel.
  norm.samples = {0.0f, 0.01f, 3.0f, 3.0f, std::numeric_limits<float>::quiet_NaN(), 3.0f, 0.0f};
  norm.original_positions = {10, 11, 12, 13, 14, 15, 16};

  TokenizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -2.0f;
  cfg.fine_max = 2.0f;
  cfg.fine_range = 0.4f;
  cfg.qbits = 4;

  auto plain = make_tokenizer(cfg)->quantize(norm);
  CHECK(plain.tokens.size() == 7);

  cfg.collapse_repeats = true;
  auto out = make_tokenizer(cfg)->quantize(norm);
  const std::int16_t sentinel = std::numeric_limits<std::int16_t>::min();
  REQUIRE(out.tokens.size() == 5);
  CHECK(out.tokens[0] == 3);
  CHECK(out.tokens[1] == 15);
  CHECK(out.tokens[2] == sentinel);
  CHECK(out.tokens[3] == 15);  // the sentinel resets the comparison
  CHECK(out.tokens[4] == 3);
  CHECK(out.original_positions == std::vector<std::uint32_t>{10, 12, 14, 15, 16});
}
