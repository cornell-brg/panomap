// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "signal/alignment_quantizers/alignment_quantizer.hpp"

namespace piru::signal {

class PassthroughAlignmentQuantizer : public AlignmentQuantizer {
public:
    explicit PassthroughAlignmentQuantizer(AlignmentQuantizerConfig config);

    AlignmentQuantizedSignal quantize(const NormalizedSignal& signal,
                                      const EventSeries* events) const override;
    const AlignmentQuantizerConfig& config() const override;
    std::string name() const override;
    float scale() const override;
    float offset() const override;

private:
    AlignmentQuantizerConfig config_;
};

}  // namespace piru::signal
