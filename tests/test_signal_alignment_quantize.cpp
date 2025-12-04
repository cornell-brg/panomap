// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"

#include <doctest/doctest.h>
#include <limits>

using namespace piru::signal;

TEST_CASE("Int16 alignment quantizer rounds and clamps") {
    NormalizedSignal norm;
    norm.samples = {1.4f, -40000.0f, 40000.0f};

    AlignmentQuantizerConfig cfg;
    cfg.backend = "int16";

    auto quantizer = make_alignment_quantizer(cfg);
    auto out = quantizer->quantize(norm, nullptr);

    CHECK(out.kind == AlignmentQuantizationKind::kInt16);
    auto* vec = std::get_if<std::vector<std::int16_t>>(&out.data);
    REQUIRE(vec != nullptr);
    REQUIRE(vec->size() == 3);
    CHECK((*vec)[0] == 1);
    CHECK((*vec)[1] == std::numeric_limits<std::int16_t>::min());
    CHECK((*vec)[2] == std::numeric_limits<std::int16_t>::max());
}

TEST_CASE("Int8 alignment quantizer outputs int8 payload") {
    NormalizedSignal norm;
    norm.samples = {5.0f, -2.4f};

    AlignmentQuantizerConfig cfg;
    cfg.backend = "int8";

    auto quantizer = make_alignment_quantizer(cfg);
    auto out = quantizer->quantize(norm, nullptr);

    CHECK(out.kind == AlignmentQuantizationKind::kInt8);
    auto* vec = std::get_if<std::vector<std::int8_t>>(&out.data);
    REQUIRE(vec != nullptr);
    REQUIRE(vec->size() == 2);
    CHECK((*vec)[0] == 5);
    CHECK((*vec)[1] == -2);
}
