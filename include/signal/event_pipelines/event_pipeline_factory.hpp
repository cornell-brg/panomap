// SPDX-License-Identifier: MIT
// Factory helpers for event pipeline construction.

#pragma once

#include "signal/event_pipeline.hpp"

namespace piru::signal {

EventPipelinePtr make_event_pipeline(const EventPipelineConfig& config);

}  // namespace piru::signal
