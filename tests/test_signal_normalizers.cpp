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

TEST_CASE("Z-score normalizer with outlier clipping") {
    piru::io::RawRead read;
    // Pattern: one large negative outlier, rest around 0
    // Values: [-1000, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    //
    // For this pattern (one outlier, n-1 zeros, n=10):
    //  - The outlier's z-score = -sqrt(n-1) ≈ -3
    //  - The others are small positive z (~ +0.33)
    // So clipping at [-2, 2] will:
    //  - Clamp the outlier to -2
    //  - Leave the rest unclipped but “near” 0.
    read.raw_signal = {-1000, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    read.range = 1.0f;
    read.digitisation = 1.0f;
    read.offset = 0.0f;

    SignalNormalizerConfig config;
    config.backend = "zscore";
    config.clip_outliers = true;
    config.clip_min = -2.0f;
    config.clip_max =  2.0f;

    auto normalizer = make_signal_normalizer(config);
    auto normalized = normalizer->normalize(read, nullptr);

    REQUIRE(normalized.samples.size() == 10);

    // Extreme outlier should be clamped to -2
    CHECK(normalized.samples[0] == doctest::Approx(config.clip_min));

    // Non-outlier points should be small z-scores (≈ +0.33), so near 0
    for (std::size_t i = 1; i < normalized.samples.size(); ++i) {
        CHECK(normalized.samples[i] == doctest::Approx(0.0f).epsilon(0.5f));
    }
}

TEST_CASE("Median-MAD normalizer with outlier clipping") {
    piru::io::RawRead read;
    // Pattern designed so that:
    //  - Median = 0
    //  - MAD = 1 (non-zero, so robust z-scores are well-defined)
    //
    // raw = {-50, -1, -1, 0, 0, 0, 0, 1, 1, 50}
    //
    // abs deviations from median (0) = [50,1,1,0,0,0,0,1,1,50]
    // sorted = [0,0,0,0,1,1,1,1,50,50]
    // MAD = 1
    //
    // So robust z for ±50 is very large (|z| >> 3), and will be clipped
    // to [-3, 3].
    read.raw_signal = {-50, -1, -1, 0, 0, 0, 0, 1, 1, 50};
    read.range = 1.0f;
    read.digitisation = 1.0f;
    read.offset = 0.0f;

    SignalNormalizerConfig config;
    config.backend = "median_mad";
    config.clip_outliers = true;
    config.clip_min = -3.0f;
    config.clip_max =  3.0f;

    auto normalizer = make_signal_normalizer(config);
    auto normalized = normalizer->normalize(read, nullptr);

    REQUIRE(normalized.samples.size() == 10);

    // Extreme outliers should be clamped to [-3, 3]
    CHECK(normalized.samples.front() == doctest::Approx(config.clip_min)); // -50
    CHECK(normalized.samples.back()  == doctest::Approx(config.clip_max)); // +50

    // Values near the median (−1, 0, +1) should map to small z-scores near 0
    for (std::size_t i = 1; i + 1 < normalized.samples.size(); ++i) {
        CHECK(normalized.samples[i] == doctest::Approx(0.0f).epsilon(0.5f));
    }
}
