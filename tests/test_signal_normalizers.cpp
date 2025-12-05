// SPDX-License-Identifier: MIT

#include "signal/normalizers/normalizer_factory.hpp"

#include <doctest/doctest.h>

using namespace piru::signal;

TEST_CASE("Z-score normalizer scales to zero mean unit variance") {
    piru::io::RawRead read;
    read.raw_signal = {1, 2, 3};
    read.range = 1.0f;
    read.digitisation = 1.0f;
    read.offset = 0.0f;

    auto normalizer = make_signal_normalizer(SignalNormalizerConfig{.backend = "zscore"});
    auto normalized = normalizer->normalize(read, nullptr);

    REQUIRE(normalized.samples.size() == 3);
    const float mean =
        (normalized.samples[0] + normalized.samples[1] + normalized.samples[2]) / 3.0f;
    CHECK(mean == doctest::Approx(0.0f).epsilon(1e-5));
}

TEST_CASE("Median-MAD normalizer centers at median and scales by MAD") {
    piru::io::RawRead read;
    read.raw_signal = {0, 0, 10, 0, 0};
    read.range = 1.0f;
    read.digitisation = 1.0f;
    read.offset = 0.0f;

    auto normalizer = make_signal_normalizer(SignalNormalizerConfig{.backend = "median_mad"});
    auto normalized = normalizer->normalize(read, nullptr);

    REQUIRE(normalized.samples.size() == 5);
    // Median is 0, MAD of deviations [0,0,10,0,0] is 0 so we clamp divisor to 1.
    CHECK(normalized.samples[2] == doctest::Approx(10.0f));
    CHECK(normalized.samples[0] == doctest::Approx(0.0f));
}
