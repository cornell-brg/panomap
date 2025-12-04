// SPDX-License-Identifier: MIT
// Factory helpers for event detector construction.

#pragma once

#include "signal/event_detectors/event_detector.hpp"

namespace piru::signal {

EventDetectorPtr make_event_detector(const EventDetectorConfig& config);

}  // namespace piru::signal
