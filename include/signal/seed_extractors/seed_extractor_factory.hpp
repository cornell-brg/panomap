// SPDX-License-Identifier: MIT
// Factory helpers for seed extractors.

#pragma once

#include "signal/seed_extractors/seed_extractor.hpp"

namespace piru::signal {

SeedExtractorPtr make_seed_extractor(const SeedExtractorConfig& config);

}  // namespace piru::signal
