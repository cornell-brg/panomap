// SPDX-License-Identifier: MIT
// Interface for detecting signal events from raw reads.

#pragma once

#include <memory>
#include <string>

#include "io/reads/read_provider.hpp"
#include "signal/signal_types.hpp"

namespace piru::signal {

struct EventDetectorConfig {
    std::string backend{"scrappie"};
    int window_length1{3};
    int window_length2{6};
    float threshold1{1.4f};
    float threshold2{9.0f};
    float peak_height{0.2f};
    int trim_start{200};
    int trim_end{10};
    int varseg_chunk{100};
    float varseg_thresh{0.0f};
};

class EventDetector {
public:
    virtual ~EventDetector() = default;

    virtual EventSeries detect(const io::RawRead& read) const = 0;
    virtual const EventDetectorConfig& config() const = 0;
    virtual std::string name() const = 0;
};

using EventDetectorPtr = std::unique_ptr<EventDetector>;

}  // namespace piru::signal
