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
  std::string backend{"rawhash"};  // "rawhash" | "scrappie" | "passthrough"
  std::string pore_model;          // e.g., "r9.4", "r10.4" - used for backend-specific defaults

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

  // Chunked processing
  std::size_t chunk_size{4000};  // Samples per chunk (0 = no chunking, process all at once)
  std::size_t max_chunks{10};    // Max chunks to process (0 = unlimited, RH2 default: 10)

  // Signal trimming parameters
  int trim_start{0};
  int trim_end{0};
  int varseg_chunk{100};
  float varseg_thresh{0.0f};

  // Normalization parameters
  std::string norm_method{"zscore"};  // "zscore" | "median_mad"
  bool clip_outliers{false};
  float clip_min{-3.0f};
  float clip_max{3.0f};
};

// Normalization state accumulated across chunks (like RH2's mean_sum/std_dev_sum).
// Allows progressive normalization: chunk N uses stats from samples 0..N.
struct NormState {
  double sum{0.0};
  double sum_sq{0.0};
  std::size_t n{0};
};

class EventPipeline {
public:
  virtual ~EventPipeline() = default;

  // Process a raw read and return normalized event signal.
  // The internal ordering (detect->normalize or normalize->detect) is backend-specific.
  virtual NormalizedSignal process(const io::RawRead& read) const = 0;

  // Process one chunk of pA signal with accumulated normalization state.
  // norm_state is updated with this chunk's statistics (progressive normalization).
  // Event detection runs independently on this chunk's normalized signal.
  virtual NormalizedSignal process_chunk(const float* pA, std::size_t len,
                                         NormState& norm_state) const = 0;

  virtual const EventPipelineConfig& config() const = 0;
  virtual std::string name() const = 0;
};

using EventPipelinePtr = std::unique_ptr<EventPipeline>;

}  // namespace piru::signal
