#include <doctest/doctest.h>

#include "core/util/metrics.hpp"

TEST_CASE("realtime is monotonic") {
  const double t1 = panomap::realtime();
  const double t2 = panomap::realtime();
  CHECK(t2 >= t1);
}

TEST_CASE("cputime is non-negative") { CHECK(panomap::cputime() >= 0.0); }
