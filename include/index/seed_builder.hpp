// SPDX-License-Identifier: MIT
// Seed extraction and store population for indexing.
//
// Design: Extensible seed extraction via SeedExtractor interface.
//
// Current support: k-mer seeds (KmerSeedExtractor)
// Future support: minimizer seeds (MinimizerSeedExtractor)
//
// To add minimizer support:
//   1. Implement MinimizerSeedExtractor : public SeedExtractor
//   2. Update make_seed_extractor factory to handle "minimizer" backend
//   3. Use: make_seed_extractor({.backend = "minimizer", .window = 10})
//   4. buildSeedStore works unchanged (backend-agnostic)

#pragma once

#include <vector>

#include "index/seed_store.hpp"
#include "signal/seed_extractors/seed_extractor.hpp"
#include "signal/signal_types.hpp"

namespace piru::index {

struct SeedBuildConfig {
    // Keep only the least frequent fraction of seeds (GraphAligner-style filtering).
    // Range: [0.0, 1.0]. Default 1.0 = keep all seeds.
    // Example: 0.1 = keep only seeds with frequency <= 10th percentile.
    double keep_least_frequent_fraction{1.0};
};

// Extract seeds from fuzzy-quantized signals and populate a SeedStore.
// Computes frequency statistics and optionally filters repetitive seeds.
HashSeedStore buildSeedStore(const std::vector<piru::signal::FuzzyQuantizedSignal>& signals,
                             const piru::signal::SeedExtractor& extractor,
                             const SeedBuildConfig& config = SeedBuildConfig{});

}  // namespace piru::index
