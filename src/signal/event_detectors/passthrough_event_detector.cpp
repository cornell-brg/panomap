// SPDX-License-Identifier: MIT

#include "signal/event_detectors/passthrough_event_detector.hpp"

namespace piru::signal {

PassthroughEventDetector::PassthroughEventDetector(EventDetectorConfig config)
    : config_(std::move(config)) {}

EventSeries PassthroughEventDetector::detect(const NormalizedSignal& signal) const {
    (void)signal;
    return {};
}

const EventDetectorConfig& PassthroughEventDetector::config() const { return config_; }

std::string PassthroughEventDetector::name() const { return config_.backend; }

}  // namespace piru::signal
