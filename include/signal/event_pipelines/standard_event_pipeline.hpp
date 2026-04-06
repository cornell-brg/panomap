// SPDX-License-Identifier: MIT
// Standard event pipeline: normalize raw signal first, then detect events.

#pragma once

#include <string>

#include "signal/event_pipeline.hpp"

namespace piru::signal {

class StandardEventPipeline : public EventPipeline {
public:
  explicit StandardEventPipeline(EventPipelineConfig config);

  NormalizedSignal process(const io::RawRead& read) const override;
  NormalizedSignal process_chunk(const float* pA, std::size_t len,
                                 NormState& norm_state) const override;
  const EventPipelineConfig& config() const override;
  std::string name() const override;

private:
  EventPipelineConfig config_;
};

}  // namespace piru::signal
