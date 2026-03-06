// SPDX-License-Identifier: MIT
// RawHash-style event pipeline: normalize raw signal first, then detect events.

#pragma once

#include <string>

#include "signal/event_pipeline.hpp"

namespace piru::signal {

class RawHashEventPipeline : public EventPipeline {
public:
  explicit RawHashEventPipeline(EventPipelineConfig config);

  NormalizedSignal process(const io::RawRead& read) const override;
  const EventPipelineConfig& config() const override;
  std::string name() const override;

private:
  EventPipelineConfig config_;
};

}  // namespace piru::signal
