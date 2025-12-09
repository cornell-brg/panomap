// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"

#include <cmath> // Added for std::isnan
#include <doctest/doctest.h>
#include <limits>

using namespace piru::signal;

TEST_CASE("Int16 alignment quantizer rounds and clamps with dynamic scale") {
    NormalizedSignal norm;
    // Input is already clipped to [-3, 3] by the normalizer
    norm.samples = {1.4f, -4.0f, 4.0f};

    AlignmentQuantizerConfig cfg;
    cfg.backend = "int16";

    auto quantizer = make_alignment_quantizer(cfg);
    auto out = quantizer->quantize(norm, nullptr);

    CHECK(out.kind == AlignmentQuantizationKind::kInt16);
    auto* vec = std::get_if<std::vector<std::int16_t>>(&out.data);
    REQUIRE(vec != nullptr);
    REQUIRE(vec->size() == 3);

    // New scale is (32766 / 3.0) = 10922.0
    // 1.4 * 10922.0 = 15290.8 -> 15291
    CHECK((*vec)[0] == 15291);
    // min() is reserved for sentinel, so clamps to min()+1
    CHECK((*vec)[1] == std::numeric_limits<std::int16_t>::min() + 1);
    CHECK((*vec)[2] == std::numeric_limits<std::int16_t>::max());
}

TEST_CASE("Int8 alignment quantizer uses dynamic scale") {
    NormalizedSignal norm;
    norm.samples = {3.0f, -2.4f};

    AlignmentQuantizerConfig cfg;
    cfg.backend = "int8";

    auto quantizer = make_alignment_quantizer(cfg);
    auto out = quantizer->quantize(norm, nullptr);

    CHECK(out.kind == AlignmentQuantizationKind::kInt8);
    auto* vec = std::get_if<std::vector<std::int8_t>>(&out.data);
    REQUIRE(vec != nullptr);
    REQUIRE(vec->size() == 2);

    // New scale is (126 / 3.0) = 42.0
    // 3.0 * 42.0 = 126
    CHECK((*vec)[0] == 126);
    // -2.4 * 42.0 = -100.8 -> -101
    CHECK((*vec)[1] == -101);
}

TEST_CASE("Passthrough alignment quantizer keeps floats") {
    NormalizedSignal norm;
    norm.samples = {1.4f, -2.8f, 0.0f, std::numeric_limits<float>::quiet_NaN()};

    AlignmentQuantizerConfig cfg;
    cfg.backend = "passthrough";

    auto quantizer = make_alignment_quantizer(cfg);
    auto out = quantizer->quantize(norm, nullptr);

    CHECK(out.kind == AlignmentQuantizationKind::kFloat32);
    auto* vec = std::get_if<std::vector<float>>(&out.data);
    REQUIRE(vec != nullptr);
    REQUIRE(vec->size() == 4);
    CHECK((*vec)[0] == 1.4f);
    CHECK((*vec)[1] == -2.8f);
    CHECK((*vec)[2] == 0.0f);
    CHECK(std::isnan((*vec)[3]));
}
