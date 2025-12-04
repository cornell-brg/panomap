// SPDX-License-Identifier: MIT

#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"

#include "signal/alignment_quantizers/passthrough_alignment_quantizer.hpp"
#include "util/logging.hpp"

namespace piru::signal {

AlignmentQuantizerPtr make_alignment_quantizer(const AlignmentQuantizerConfig& config) {
    if (config.backend != "int16" && config.backend != "int8" && config.backend != "passthrough" &&
        config.backend != "float32") {
        LOG_WARN("Unknown alignment quantizer backend '" + config.backend + "', using passthrough");
    }
    return std::make_unique<PassthroughAlignmentQuantizer>(config);
}

}  // namespace piru::signal
