// SPDX-License-Identifier: MIT

#include "signal/event_detectors/passthrough_event_detector.hpp"

namespace piru::signal {

PassthroughEventDetector::PassthroughEventDetector(EventDetectorConfig config)
    : config_(std::move(config)) {}

EventSeries PassthroughEventDetector::detect(const io::RawRead& read) const {
    EventSeries series;
    series.sampling_rate_hz = read.sampling_rate_hz;

    // Convert each sample to an event (one sample per event)
    // Convert ADC values to picoamps
    const float digitisation = read.digitisation;
    const float offset = read.offset;
    const float range = read.range;

    const std::size_t num_samples = read.raw_signal.size();
    series.events.reserve(num_samples);
    for (std::size_t i = 0; i < num_samples; ++i) {
        const float pA = ((read.raw_signal[i] + offset) * range) / digitisation;
        series.events.push_back(SignalEvent{
            .start = i,
            .length = 1,
            .mean = pA,
            .stdv = 0.0f
        });
    }

    return series;
}

const EventDetectorConfig& PassthroughEventDetector::config() const { return config_; }

std::string PassthroughEventDetector::name() const { return config_.backend; }

}  // namespace piru::signal
