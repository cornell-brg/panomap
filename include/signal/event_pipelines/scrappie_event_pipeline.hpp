// SPDX-License-Identifier: MIT
// Scrappie-style event pipeline: detect events on raw signal, then normalize.

#pragma once

#include <string>

#include "signal/event_pipeline.hpp"

namespace piru::signal {

class ScrappieEventPipeline : public EventPipeline {
public:
  explicit ScrappieEventPipeline(EventPipelineConfig config);

  NormalizedSignal process(const io::RawRead& read) const override;
  const EventPipelineConfig& config() const override;
  std::string name() const override;

private:
  EventPipelineConfig config_;
};

}  // namespace piru::signal
