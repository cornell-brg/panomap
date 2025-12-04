// SPDX-License-Identifier: MIT

#include "signal/event_detectors/event_detector_factory.hpp"

#include <doctest/doctest.h>

using namespace piru::signal;

TEST_CASE("Scrappie detector segments simple step signal") {
    NormalizedSignal norm;
    norm.samples.assign(10, 0.0f);
    norm.samples.insert(norm.samples.end(), 10, 10.0f);

    EventDetectorConfig cfg;
    cfg.backend = "scrappie";
    cfg.trim_start = 0;
    cfg.trim_end = 0;
    cfg.varseg_chunk = 100;  // larger than signal to avoid trimming
    cfg.varseg_thresh = 0.0f;
    cfg.threshold1 = 0.5f;
    cfg.threshold2 = 1.0f;

    auto detector = make_event_detector(cfg);
    auto events = detector->detect(norm);

    REQUIRE(events.events.size() == 2);
    CHECK(events.events[0].start == 0);
    CHECK(events.events[0].length == 10);
    CHECK(events.events[0].mean == doctest::Approx(0.0f));
    CHECK(events.events[1].start == 10);
    CHECK(events.events[1].length == 10);
    CHECK(events.events[1].mean == doctest::Approx(10.0f));
}

TEST_CASE("Scrappie detector segments multiple transitions") {
    NormalizedSignal norm;
    norm.samples.insert(norm.samples.end(), 5, 0.0f);
    norm.samples.insert(norm.samples.end(), 5, 10.0f);
    norm.samples.insert(norm.samples.end(), 5, -5.0f);

    EventDetectorConfig cfg;
    cfg.backend = "scrappie";
    cfg.trim_start = 0;
    cfg.trim_end = 0;
    cfg.varseg_chunk = 100;
    cfg.varseg_thresh = 0.0f;
    cfg.threshold1 = 0.5f;
    cfg.threshold2 = 1.0f;

    auto detector = make_event_detector(cfg);
    auto events = detector->detect(norm);

    REQUIRE(events.events.size() == 3);
    CHECK(events.events[0].mean == doctest::Approx(0.0f));
    CHECK(events.events[1].mean == doctest::Approx(10.0f));
    CHECK(events.events[2].mean == doctest::Approx(-5.0f));
    CHECK(events.events[0].length == 5);
    CHECK(events.events[1].length == 5);
    CHECK(events.events[2].length == 5);
}

TEST_CASE("Scrappie detector handles empty input") {
    NormalizedSignal norm;
    EventDetectorConfig cfg;
    cfg.backend = "scrappie";
    cfg.trim_start = 0;
    cfg.trim_end = 0;

    auto detector = make_event_detector(cfg);
    auto events = detector->detect(norm);
    CHECK(events.events.empty());
}
