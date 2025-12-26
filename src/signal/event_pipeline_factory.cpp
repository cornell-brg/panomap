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
    LOG_INFO("Event pipeline: user overrides applied (w1=" + std::to_string(config.window_length1) +
             ", w2=" + std::to_string(config.window_length2) +
             ", t1=" + std::to_string(config.threshold1) +
             ", t2=" + std::to_string(config.threshold2) +
             ", ph=" + std::to_string(config.peak_height) + ")");
  }
}

// Apply RawHash-specific defaults based on pore model (R9 vs R10)
EventPipelineConfig apply_rawhash_defaults(EventPipelineConfig config) {
  const bool is_r10 = (config.pore_model.find("r10") != std::string::npos ||
                       config.pore_model.find("R10") != std::string::npos);

  if (is_r10) {
    // R10 parameters (from RawHash main.cpp preset)
    config.window_length1 = 3;
    config.window_length2 = 6;
    config.threshold1 = 6.5f;
    config.threshold2 = 4.0f;
    config.peak_height = 0.2f;
    LOG_INFO("RawHash event pipeline: using R10 defaults (w1=3, w2=6, t1=6.5, t2=4.0, ph=0.2)");
  } else {
    // R9 parameters (RawHash defaults from roptions.c)
    config.window_length1 = 3;
    config.window_length2 = 9;
    config.threshold1 = 4.0f;
    config.threshold2 = 3.5f;
    config.peak_height = 0.4f;
    LOG_INFO("RawHash event pipeline: using R9 defaults (w1=3, w2=9, t1=4.0, t2=3.5, ph=0.4)");
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
  if (config.backend == "scrappie" || config.backend.empty()) {
    LOG_INFO("Event pipeline: scrappie (detect events on raw signal, then normalize)");
    auto scrappie_config = apply_scrappie_defaults(config);
    return std::make_unique<ScrappieEventPipeline>(scrappie_config);
  }
  if (config.backend == "rawhash") {
    LOG_INFO("Event pipeline: rawhash (normalize raw signal, then detect events)");
    auto rawhash_config = apply_rawhash_defaults(config);
    return std::make_unique<RawHashEventPipeline>(rawhash_config);
  }
  if (config.backend == "passthrough") {
    LOG_INFO("Event pipeline: passthrough (normalize only, no event detection)");
    return std::make_unique<PassthroughEventPipeline>(config);
  }

  LOG_WARN("Unknown event pipeline backend '" + config.backend + "', using scrappie");
  auto scrappie_config = apply_scrappie_defaults(config);
  return std::make_unique<ScrappieEventPipeline>(scrappie_config);
}

}  // namespace piru::signal
