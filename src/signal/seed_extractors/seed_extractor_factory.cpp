// SPDX-License-Identifier: MIT

#include "signal/seed_extractors/seed_extractor_factory.hpp"

#include "signal/seed_extractors/kmer_seed_extractor.hpp"
#include "signal/seed_extractors/minimizer_seed_extractor.hpp"
#include "util/logging.hpp"

namespace piru::signal {

SeedExtractorPtr make_seed_extractor(const SeedExtractorConfig& config) {
    if (config.backend == "minimizer") {
        return std::make_unique<MinimizerSeedExtractor>(config);
    }
    if (config.backend != "kmer") {
        LOG_WARN("Unknown seed extractor backend '" + config.backend + "', using kmer");
    }
    return std::make_unique<KmerSeedExtractor>(config);
}

}  // namespace piru::signal
