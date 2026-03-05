// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/event_pipeline_factory.hpp"

#include "signal/event_pipelines/passthrough_event_pipeline.hpp"
#include "signal/event_pipelines/rawhash_event_pipeline.hpp"
#include "signal/event_pipelines/scrappie_event_pipeline.hpp"
#include "util/logging.hpp"

namespace piru::signal {

namespace {

// Apply user overrides to config (if any were set via CLI)
void apply_user_overrides(EventPipelineConfig& config) {
    bool has_overrides = false;
    if (config.override_window_length1) {
        config.window_length1 = *config.override_window_length1;
        has_overrides = true;
    }
    if (config.override_window_length2) {
        config.window_length2 = *config.override_window_length2;
        has_overrides = true;
    }
    if (config.override_threshold1) {
        config.threshold1 = *config.override_threshold1;
        has_overrides = true;
    }
    if (config.override_threshold2) {
        config.threshold2 = *config.override_threshold2;
        has_overrides = true;
    }
    if (config.override_peak_height) {
        config.peak_height = *config.override_peak_height;
        has_overrides = true;
    }
    if (has_overrides) {
        LOG_DEBUG(
            "Event pipeline: user overrides applied (w1=" + std::to_string(config.window_length1) +
            ", w2=" + std::to_string(config.window_length2) + ", t1=" +
            std::to_string(config.threshold1) + ", t2=" + std::to_string(config.threshold2) +
            ", ph=" + std::to_string(config.peak_height) + ")");
    }
}

// Apply RawHash-specific defaults based on pore model (R9 vs R10)
EventPipelineConfig apply_rawhash_defaults(EventPipelineConfig config) {
    const bool is_r10 = (config.pore_model.find("r10") != std::string::npos ||
                         config.pore_model.find("R10") != std::string::npos);

    if (is_r10) {
        // R10 parameters (tuned for piru, see DEV025)
        config.window_length1 = 4;
        config.window_length2 = 10;
        config.threshold1 = 4.0f;
        config.threshold2 = 3.0f;
        config.peak_height = 0.4f;
        LOG_DEBUG(
            "RawHash event pipeline: using R10 defaults (w1=4, w2=10, t1=4.0, t2=3.0, ph=0.4)");
        // Original RawHash R10 defaults (from RawHash main.cpp preset):
        // w1=3, w2=6, t1=6.5, t2=4.0, ph=0.2
    } else {
        // R9 parameters (tuned for piru, see DEV025)
        config.window_length1 = 5;
        config.window_length2 = 8;
        config.threshold1 = 2.1f;
        config.threshold2 = 1.2f;
        config.peak_height = 0.4f;
        LOG_DEBUG("RawHash event pipeline: using R9 defaults (w1=5, w2=8, t1=2.1, t2=1.2, ph=0.4)");
        // Original RawHash R9 defaults (from RawHash roptions.c):
        // w1=3, w2=9, t1=4.0, t2=3.5, ph=0.4
    }

    // Apply user overrides after backend defaults
    apply_user_overrides(config);

    return config;
}

// Apply Scrappie defaults (and user overrides)
EventPipelineConfig apply_scrappie_defaults(EventPipelineConfig config) {
    // Scrappie uses the struct defaults, just apply overrides
    apply_user_overrides(config);
    return config;
}

}  // namespace

EventPipelinePtr make_event_pipeline(const EventPipelineConfig& config) {
    if (config.backend == "rawhash" || config.backend.empty()) {
        LOG_DEBUG("Event pipeline: rawhash (normalize raw signal, then detect events)");
        auto rawhash_config = apply_rawhash_defaults(config);
        return std::make_unique<RawHashEventPipeline>(rawhash_config);
    }
    if (config.backend == "scrappie") {
        LOG_DEBUG("Event pipeline: scrappie (detect events on raw signal, then normalize)");
        auto scrappie_config = apply_scrappie_defaults(config);
        return std::make_unique<ScrappieEventPipeline>(scrappie_config);
    }
    if (config.backend == "passthrough") {
        LOG_DEBUG("Event pipeline: passthrough (normalize only, no event detection)");
        return std::make_unique<PassthroughEventPipeline>(config);
    }

    LOG_WARN("Unknown event pipeline backend '" + config.backend + "', using rawhash");
    auto rawhash_config = apply_rawhash_defaults(config);
    return std::make_unique<RawHashEventPipeline>(rawhash_config);
}

}  // namespace piru::signal
