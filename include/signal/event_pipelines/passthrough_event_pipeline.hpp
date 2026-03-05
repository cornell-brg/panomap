// SPDX-License-Identifier: MIT
// Passthrough event pipeline: normalize raw samples without event detection.

#pragma once

#include <string>

#include "signal/event_pipeline.hpp"

namespace piru::signal {

class PassthroughEventPipeline : public EventPipeline {
public:
    explicit PassthroughEventPipeline(EventPipelineConfig config);

    NormalizedSignal process(const io::RawRead& read) const override;
    const EventPipelineConfig& config() const override;
    std::string name() const override;

private:
    EventPipelineConfig config_;
};

}  // namespace piru::signal
