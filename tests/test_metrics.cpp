#include <doctest/doctest.h>

#include "util/metrics.hpp"

TEST_CASE("realtime is monotonic") {
    const double t1 = piru::realtime();
    const double t2 = piru::realtime();
    CHECK(t2 >= t1);
}

TEST_CASE("cputime is non-negative") { CHECK(piru::cputime() >= 0.0); }
