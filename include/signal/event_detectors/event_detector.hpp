// SPDX-License-Identifier: MIT
// Interface for detecting signal events from normalized reads.

#pragma once

#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct EventDetectorConfig {
    std::string backend{"scrappie"};
};

class EventDetector {
public:
    virtual ~EventDetector() = default;

    virtual EventSeries detect(const NormalizedSignal& signal) const = 0;
    virtual const EventDetectorConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using EventDetectorPtr = std::unique_ptr<EventDetector>;

}  // namespace piru::signal
