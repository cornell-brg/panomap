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

TEST_CASE("RH2 bucket quantizer gives ordered buckets that use every value") {
  TokenizerConfig cfg;
  cfg.backend = "rh2";
  cfg.fine_min = -1.6f;
  cfg.fine_max = 1.6f;
  cfg.fine_range = 0.6f;
  cfg.n_buckets = 12;

  NormalizedSignal norm;
  for (int i = 0; i <= 600; ++i) norm.samples.push_back(-3.0f + 0.01f * static_cast<float>(i));
  auto out = make_tokenizer(cfg)->quantize(norm);

  REQUIRE(out.tokens.size() == norm.samples.size());
  CHECK(out.tokens.front() == 0);
  CHECK(out.tokens.back() == 11);
  for (std::size_t i = 1; i < out.tokens.size(); ++i) CHECK(out.tokens[i] >= out.tokens[i - 1]);
  // 12 buckets: 7 in the fine range, 2 below it, 3 above it.
  CHECK(make_tokenizer(cfg)->quantize(NormalizedSignal{{-1.7f}}).tokens[0] == 1);
  CHECK(make_tokenizer(cfg)->quantize(NormalizedSignal{{-1.5f}}).tokens[0] == 2);
  CHECK(make_tokenizer(cfg)->quantize(NormalizedSignal{{1.5f}}).tokens[0] == 8);
  CHECK(make_tokenizer(cfg)->quantize(NormalizedSignal{{1.7f}}).tokens[0] == 9);
}

TEST_CASE("Token bits follow the tokenizer") {
  TokenizerConfig cfg;
  CHECK(token_bits(cfg) == 4);
  cfg.n_buckets = 24;
  CHECK(token_bits(cfg) == 5);
  cfg.backend = "landmark";
  cfg.landmark_bits = 3;
  CHECK(token_bits(cfg) == 6);
}
