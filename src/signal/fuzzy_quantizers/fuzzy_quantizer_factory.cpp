// SPDX-License-Identifier: MIT

#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"

#include "signal/fuzzy_quantizers/passthrough_fuzzy_quantizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

FuzzyQuantizerPtr make_fuzzy_quantizer(const FuzzyQuantizerConfig& config) {
    if (config.backend != "rawhash2" && config.backend != "passthrough") {
        LOG_WARN("Unknown fuzzy quantizer backend '" + config.backend + "', using passthrough");
    }
    return std::make_unique<PassthroughFuzzyQuantizer>(config);
}

}  // namespace piru::signal
