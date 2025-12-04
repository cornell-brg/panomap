// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/kmer_seed_extractor.hpp"

namespace piru::signal {

KmerSeedExtractor::KmerSeedExtractor(SeedExtractorConfig config) : config_(std::move(config)) {}

SeedBuffer KmerSeedExtractor::extract(const FuzzyQuantizedSignal& signal,
                                      const EventSeries* events) const {
    (void)events;
    SeedBuffer buffer;
    buffer.seeds.reserve(signal.tokens.size());
    for (std::size_t i = 0; i < signal.tokens.size(); ++i) {
        const auto token = signal.tokens[i];
        buffer.seeds.push_back(
            Seed{.hash = static_cast<std::uint64_t>(token), .position = i, .span = 1});
    }
    return buffer;
}

const SeedExtractorConfig& KmerSeedExtractor::config() const { return config_; }

std::string KmerSeedExtractor::name() const { return config_.backend; }

}  // namespace piru::signal
