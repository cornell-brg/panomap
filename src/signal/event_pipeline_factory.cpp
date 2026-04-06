// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/event_pipeline_factory.hpp"

#include "signal/event_pipelines/standard_event_pipeline.hpp"
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
        ", w2=" + std::to_string(config.window_length2) +
        ", t1=" + std::to_string(config.threshold1) + ", t2=" + std::to_string(config.threshold2) +
        ", ph=" + std::to_string(config.peak_height) + ")");
  }
}

// Apply pore model-specific defaults (R9 vs R10)
EventPipelineConfig apply_defaults(EventPipelineConfig config) {
  const bool is_r10 = (config.pore_model.find("r10") != std::string::npos ||
                       config.pore_model.find("R10") != std::string::npos);

  if (is_r10) {
    config.window_length1 = 4;
    config.window_length2 = 10;
    config.threshold1 = 4.0f;
    config.threshold2 = 3.0f;
    config.peak_height = 0.4f;
    LOG_DEBUG("Event pipeline: using R10 defaults (w1=4, w2=10, t1=4.0, t2=3.0, ph=0.4)");
  } else {
    config.window_length1 = 3;
    config.window_length2 = 9;
    config.threshold1 = 4.0f;
    config.threshold2 = 3.5f;
    config.peak_height = 0.4f;
    LOG_DEBUG("Event pipeline: using R9 defaults (w1=3, w2=9, t1=4.0, t2=3.5, ph=0.4)");
  }

  // Apply user overrides after defaults
  apply_user_overrides(config);

  return config;
}

}  // namespace

EventPipelinePtr make_event_pipeline(const EventPipelineConfig& config) {
  if (config.backend != "standard" && !config.backend.empty()) {
    LOG_WARN("Unknown event pipeline backend '" + config.backend + "', using standard");
  }

  auto resolved = apply_defaults(config);
  return std::make_unique<StandardEventPipeline>(resolved);
}

}  // namespace piru::signal
