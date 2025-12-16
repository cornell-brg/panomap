// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "signal/seed_extractors/seed_extractor.hpp"

namespace piru::signal {

class KmerSeedExtractor : public SeedExtractor {
public:
    explicit KmerSeedExtractor(SeedExtractorConfig config);

    SeedBuffer extract(const FuzzyQuantizedSignal& signal) const override;
    const SeedExtractorConfig& config() const override;
    std::string name() const override;

private:
    SeedExtractorConfig config_;
};

}  // namespace piru::signal
