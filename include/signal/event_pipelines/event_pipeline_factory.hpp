// SPDX-License-Identifier: MIT
// Factory helpers for event pipeline construction.

#pragma once

#include "signal/event_pipeline.hpp"

namespace panomap::signal {

EventPipelinePtr make_event_pipeline(const EventPipelineConfig& config);

}  // namespace panomap::signal
