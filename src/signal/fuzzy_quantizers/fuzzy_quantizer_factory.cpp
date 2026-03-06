// SPDX-License-Identifier: MIT

#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"

#include "signal/fuzzy_quantizers/passthrough_fuzzy_quantizer.hpp"
#include "signal/fuzzy_quantizers/piru_fuzzy_quantizer.hpp"
#include "signal/fuzzy_quantizers/rh2_fuzzy_quantizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

namespace {

// Apply chemistry-specific defaults for piru backend (see DEV025)
FuzzyQuantizerConfig apply_piru_defaults(FuzzyQuantizerConfig config) {
  const bool is_r10 = (config.pore_model.find("r10") != std::string::npos ||
                       config.pore_model.find("R10") != std::string::npos);

  if (is_r10) {
    // R10 parameters (tuned for piru, see DEV025)
    config.fine_range = 0.8f;
    config.n_bins = 10;
    LOG_INFO("Fuzzy quantizer: using R10 defaults (fine_range=0.8, n_bins=10)");
    // Original defaults: fine_range=0.4, n_bins=16
  } else {
    // R9 parameters (tuned for piru, see DEV025)
    config.fine_range = 0.9f;
    config.n_bins = 10;
    LOG_INFO("Fuzzy quantizer: using R9 defaults (fine_range=0.9, n_bins=10)");
    // Original defaults: fine_range=0.4, n_bins=16
  }

  return config;
}

}  // namespace

FuzzyQuantizerPtr make_fuzzy_quantizer(const FuzzyQuantizerConfig& config) {
  if (config.backend == "rawhash2" || config.backend == "rh2") {
    return std::make_unique<Rh2FuzzyQuantizer>(config);
  }
  if (config.backend == "piru" || config.backend.empty()) {
    auto piru_config = apply_piru_defaults(config);
    return std::make_unique<PiruFuzzyQuantizer>(piru_config);
  }
  if (config.backend != "passthrough") {
    LOG_WARN("Unknown fuzzy quantizer backend '" + config.backend + "', using passthrough");
  }
  return std::make_unique<PassthroughFuzzyQuantizer>(config);
}

}  // namespace piru::signal
