// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "signal/tokenizers/tokenizer_factory.hpp"

using namespace piru::signal;

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
