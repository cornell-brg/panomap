// SPDX-License-Identifier: MIT
// Events from given borders: one event per stretch between borders.

#include <doctest/doctest.h>

#include <vector>

#include "signal/event_pipelines/standard_event_pipeline.hpp"

using panomap::signal::events_from_borders;
using panomap::signal::NormState;

TEST_CASE("events_from_borders: one event per step, in order") {
  // Three flat steps of 20 samples: 80, 100, 120 pA.
  std::vector<float> pA;
  for (float level : {80.0f, 100.0f, 120.0f}) pA.insert(pA.end(), 20, level);
  NormState state;
  auto events = events_from_borders(pA.data(), pA.size(), state, {0, 20, 40, 60});
  REQUIRE(events.samples.size() == 3);
  CHECK(events.samples[0] < events.samples[1]);
  CHECK(events.samples[1] < events.samples[2]);
  CHECK(events.samples[1] == doctest::Approx(0.0f).epsilon(0.01));
}

TEST_CASE("events_from_borders: empty and over-long stretches give no event") {
  std::vector<float> pA(1200, 100.0f);
  for (std::size_t i = 0; i < pA.size(); ++i) pA[i] += static_cast<float>(i % 7);
  NormState state;
  // [10,10) is empty, [10,30) is kept, [30,700) is 670 samples and is dropped.
  auto events = events_from_borders(pA.data(), pA.size(), state, {10, 10, 30, 700});
  CHECK(events.samples.size() == 1);
}
