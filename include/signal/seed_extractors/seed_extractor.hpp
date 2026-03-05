// SPDX-License-Identifier: MIT
// Interface for extracting seeds from fuzzy-quantized signals.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct SeedExtractorConfig {
    std::string backend{"kmer"};
    std::size_t k{10};
    std::size_t stride{1};
    std::uint32_t qbits{4};  // Bits per token in fuzzy quantization.
    std::size_t window{0};   // Only used for minimizer-style backends.
    std::string params;
    bool operator==(const SeedExtractorConfig& other) const = default;
};

class SeedExtractor {
public:
    virtual ~SeedExtractor() = default;

    virtual SeedBuffer extract(const FuzzyQuantizedSignal& signal) const = 0;
    virtual const SeedExtractorConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using SeedExtractorPtr = std::unique_ptr<SeedExtractor>;

}  // namespace piru::signal
