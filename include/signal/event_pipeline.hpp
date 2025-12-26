// SPDX-License-Identifier: MIT
// Unified interface for event detection and normalization.
//
// This interface consolidates event detection and normalization into a single
// pipeline stage. Different backends can implement different internal orderings:
// - Scrappie-style: detect events on raw signal, then normalize
// - RawHash-style: normalize raw signal first, then detect events

#pragma once

#include <memory>
#include <optional>
#include <string>

#include "io/reads/read_provider.hpp"
#include "signal/signal_types.hpp"

namespace piru::signal {

struct EventPipelineConfig {
  std::string backend{"scrappie"};  // "scrappie" | "rawhash" | "passthrough"
  std::string pore_model;           // e.g., "r9.4", "r10.4" - used for backend-specific defaults

  // Event detection parameters (t-stat peak detection)
  // Note: For rawhash backend, these are auto-set based on pore_model in factory
  int window_length1{3};
  int window_length2{6};
  float threshold1{1.4f};
  float threshold2{9.0f};
  float peak_height{0.2f};

  // User overrides (if set, take precedence over backend defaults)
  std::optional<int> override_window_length1;
  std::optional<int> override_window_length2;
  std::optional<float> override_threshold1;
  std::optional<float> override_threshold2;
  std::optional<float> override_peak_height;

  // Signal trimming parameters
  int trim_start{200};
  int trim_end{10};
  int varseg_chunk{100};
  float varseg_thresh{0.0f};

  // Normalization parameters
  std::string norm_method{"zscore"};  // "zscore" | "median_mad"
  bool clip_outliers{false};
  float clip_min{-3.0f};
  float clip_max{3.0f};
};

class EventPipeline {
 public:
  virtual ~EventPipeline() = default;

  // Process a raw read and return normalized event signal.
  // The internal ordering (detect→normalize or normalize→detect) is backend-specific.
  virtual NormalizedSignal process(const io::RawRead& read) const = 0;

  virtual const EventPipelineConfig& config() const = 0;
  virtual std::string name() const = 0;
};

using EventPipelinePtr = std::unique_ptr<EventPipeline>;

}  // namespace piru::signal
