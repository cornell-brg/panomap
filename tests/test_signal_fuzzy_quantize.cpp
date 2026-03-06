// SPDX-License-Identifier: MIT

#include <doctest/doctest.h>

#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"

using namespace piru::signal;

TEST_CASE("RH2 quantizer bins within fine range") {
  NormalizedSignal norm;
  norm.samples = {0.0f};

  FuzzyQuantizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -2.0f;
  cfg.fine_max = 2.0f;
  cfg.fine_range = 0.4f;
  cfg.qbits = 4;

  auto quantizer = make_fuzzy_quantizer(cfg);
  auto out = quantizer->quantize(norm);

  REQUIRE(out.tokens.size() == 1);
  CHECK(out.tokens[0] == 3);  // rawhash2 formula -> bucket ~3 for 0.0
}

TEST_CASE("RH2 quantizer clamps and maps extremes") {
  NormalizedSignal norm;
  norm.samples = {-3.0f, 3.0f};

  FuzzyQuantizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -2.0f;
  cfg.fine_max = 2.0f;
  cfg.fine_range = 0.4f;
  cfg.qbits = 4;

  auto quantizer = make_fuzzy_quantizer(cfg);
  auto out = quantizer->quantize(norm);

  REQUIRE(out.tokens.size() == 2);
  CHECK(out.tokens[0] == 6);   // lower extreme maps to coarse bucket
  CHECK(out.tokens[1] == 15);  // upper extreme maps to top bucket
}
