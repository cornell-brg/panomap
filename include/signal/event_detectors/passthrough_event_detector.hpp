// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include "signal/event_detectors/event_detector.hpp"

namespace piru::signal {

class PassthroughEventDetector : public EventDetector {
public:
    explicit PassthroughEventDetector(EventDetectorConfig config);

    EventSeries detect(const NormalizedSignal& signal) const override;
    const EventDetectorConfig& config() const override;
    std::string name() const override;

private:
    EventDetectorConfig config_;
};

}  // namespace piru::signal
