// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/passthrough_event_pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace piru::signal {

namespace {

constexpr float kMadScale = 1.4826f;

float median_inplace(std::vector<float>& vec) {
  const size_t n = vec.size();
  if (n == 0) return 0.0f;
  const size_t half = n / 2;
  std::nth_element(vec.begin(), vec.begin() + half, vec.end());
  float med = vec[half];
  if ((n % 2) == 0) {
    std::nth_element(vec.begin(), vec.begin() + half - 1, vec.end());
    med = 0.5f * (med + vec[half - 1]);
  }
  return med;
}

}  // namespace

PassthroughEventPipeline::PassthroughEventPipeline(EventPipelineConfig config)
    : config_(std::move(config)) {}

NormalizedSignal PassthroughEventPipeline::process(const io::RawRead& read) const {
  NormalizedSignal result;
  result.sampling_rate_hz = read.sampling_rate_hz;

  if (read.raw_signal.empty()) {
    return result;
  }

  // Step 1: Convert raw ADC values to picoamps
  const float raw_unit = (read.digitisation == 0.0f) ? 1.0f : (read.range / read.digitisation);
  std::vector<float> values;
  values.reserve(read.raw_signal.size());
  for (auto value : read.raw_signal) {
    values.push_back((static_cast<float>(value) + read.offset) * raw_unit);
  }

  // Step 2: Normalize (no event detection, just normalize raw samples)
  if (config_.norm_method == "median_mad") {
    std::vector<float> scratch = values;
    const float med = median_inplace(scratch);
    for (auto& v : scratch) {
      v = std::abs(v - med);
    }
    float mad_val = median_inplace(scratch) * kMadScale;
    if (mad_val == 0.0f) {
      mad_val = 1.0f;
    }

    result.samples.reserve(values.size());
    for (const auto v : values) {
      float norm = (v - med) / mad_val;
      if (config_.clip_outliers) {
        norm = std::clamp(norm, config_.clip_min, config_.clip_max);
      }
      result.samples.push_back(norm);
    }
  } else {
    // Default: zscore
    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    const double mean = sum / static_cast<double>(values.size());

    double var = 0.0;
    for (const auto v : values) {
      const double diff = static_cast<double>(v) - mean;
      var += diff * diff;
    }
    var /= static_cast<double>(values.size());
    const double stddev = (var > 0.0) ? std::sqrt(var) : 1.0;

    result.samples.reserve(values.size());
    for (const auto v : values) {
      float z = static_cast<float>((static_cast<double>(v) - mean) / stddev);
      if (config_.clip_outliers) {
        z = std::clamp(z, config_.clip_min, config_.clip_max);
      }
      result.samples.push_back(z);
    }
  }

  return result;
}

const EventPipelineConfig& PassthroughEventPipeline::config() const { return config_; }

std::string PassthroughEventPipeline::name() const { return "passthrough"; }

}  // namespace piru::signal
