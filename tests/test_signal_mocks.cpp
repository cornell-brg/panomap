// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/normalizers/normalizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"

#include <doctest/doctest.h>

using namespace piru::signal;

TEST_CASE("Identity normalizer copies raw signal to floats") {
    piru::io::RawRead read;
    read.raw_signal = {1, -2, 3};
    read.sampling_rate_hz = 4000.0f;
    read.range = 1.0f;
    read.digitisation = 1.0f;

    auto normalizer = make_signal_normalizer(SignalNormalizerConfig{.backend = "identity"});
    auto normalized = normalizer->normalize(read);

    REQUIRE(normalized.samples.size() == 3);
    CHECK(normalized.samples[0] == doctest::Approx(1.0f));
    CHECK(normalized.samples[1] == doctest::Approx(-2.0f));
    CHECK(normalized.samples[2] == doctest::Approx(3.0f));
    CHECK(normalized.sampling_rate_hz == doctest::Approx(4000.0f));
}

TEST_CASE("Passthrough fuzzy quantizer truncates to int16 tokens") {
    NormalizedSignal norm;
    norm.samples = {1.7f, -2.5f, 0.2f};

    auto quantizer = make_fuzzy_quantizer(FuzzyQuantizerConfig{.backend = "passthrough"});
    auto quantized = quantizer->quantize(norm, nullptr);

    REQUIRE(quantized.tokens.size() == 3);
    // static_cast<int16_t> truncates toward zero
    CHECK(quantized.tokens[0] == static_cast<std::int16_t>(1));
    CHECK(quantized.tokens[1] == static_cast<std::int16_t>(-2));
    CHECK(quantized.tokens[2] == static_cast<std::int16_t>(0));
}

TEST_CASE("Passthrough alignment quantizer emits int16 payload") {
    NormalizedSignal norm;
    norm.samples = {5.9f, -1.1f};

    auto quantizer = make_alignment_quantizer(AlignmentQuantizerConfig{.backend = "passthrough"});
    auto quantized = quantizer->quantize(norm, nullptr);

    CHECK(quantized.kind == AlignmentQuantizationKind::kInt16);
    auto* vec_ptr = std::get_if<std::vector<std::int16_t>>(&quantized.data);
    REQUIRE(vec_ptr != nullptr);
    REQUIRE(vec_ptr->size() == 2);
    CHECK((*vec_ptr)[0] == static_cast<std::int16_t>(5));
    CHECK((*vec_ptr)[1] == static_cast<std::int16_t>(-1));
}
