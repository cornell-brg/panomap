// SPDX-License-Identifier: MIT

#include "signal/event_detectors/event_detector_factory.hpp"

#include "signal/event_detectors/passthrough_event_detector.hpp"
#include "util/logging.hpp"

namespace piru::signal {

EventDetectorPtr make_event_detector(const EventDetectorConfig& config) {
    if (config.backend != "scrappie" && config.backend != "passthrough") {
        LOG_WARN("Unknown event detector backend '" + config.backend + "', using passthrough");
    }
    return std::make_unique<PassthroughEventDetector>(config);
}

}  // namespace piru::signal
