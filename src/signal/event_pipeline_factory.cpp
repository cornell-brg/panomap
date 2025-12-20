// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/event_pipeline_factory.hpp"

#include "signal/event_pipelines/passthrough_event_pipeline.hpp"
#include "signal/event_pipelines/rawhash_event_pipeline.hpp"
#include "signal/event_pipelines/scrappie_event_pipeline.hpp"
#include "util/logging.hpp"

namespace piru::signal {

namespace {

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
    LOG_INFO("RawHash event pipeline: using R10 parameters (w1=3, w2=6, t1=6.5, t2=4.0, ph=0.2)");
  } else {
    // R9 parameters (RawHash defaults from roptions.c)
    config.window_length1 = 3;
    config.window_length2 = 9;
    config.threshold1 = 4.0f;
    config.threshold2 = 3.5f;
    config.peak_height = 0.4f;
    LOG_INFO("RawHash event pipeline: using R9 parameters (w1=3, w2=9, t1=4.0, t2=3.5, ph=0.4)");
  }

  return config;
}

}  // namespace

EventPipelinePtr make_event_pipeline(const EventPipelineConfig& config) {
  if (config.backend == "scrappie" || config.backend.empty()) {
    LOG_INFO("Event pipeline: scrappie (detect events on raw signal, then normalize)");
    return std::make_unique<ScrappieEventPipeline>(config);
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
  return std::make_unique<ScrappieEventPipeline>(config);
}

}  // namespace piru::signal
