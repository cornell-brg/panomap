// SPDX-License-Identifier: MIT
// Standard event pipeline: normalize raw signal first, then detect events.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "signal/event_pipeline.hpp"

namespace panomap::signal {

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

// Events from given borders instead of t-test peaks. `borders` are sample
// indices into this chunk, ascending. Each event is the IQR-filtered mean of the
// samples between two neighboring borders, on the same progressive normalization
// and |z| < 3 filter as process_chunk. Stretches with no kept sample, or of 500
// samples or more, give no event.
NormalizedSignal events_from_borders(const float* pA, std::size_t len, NormState& norm_state,
                                     const std::vector<std::uint32_t>& borders);

}  // namespace panomap::signal
