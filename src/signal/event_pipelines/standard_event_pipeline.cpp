// SPDX-License-Identifier: MIT

#include "signal/event_pipelines/standard_event_pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace piru::signal {

namespace {

// Standard event pipeline: normalize raw signal first, then detect events.
// 1. Zscore normalize raw samples, filter |z| >= 3
// 2. Run t-stat event detection on normalized signal
// 3. Use IQR-filtered mean for each event segment
// 4. Skip events > 500 samples

void compute_sum_sumsq(const std::vector<float>& signal, std::vector<float>& sum,
                       std::vector<float>& sumsq) {
  // Uses float for prefix sums
  const int n = static_cast<int>(signal.size());
  sum.assign(n + 1, 0.0f);
  sumsq.assign(n + 1, 0.0f);
  for (int i = 0; i < n; ++i) {
    sum[i + 1] = sum[i] + signal[i];
    sumsq[i + 1] = sumsq[i] + signal[i] * signal[i];
  }
}

std::vector<float> compute_tstat(const std::vector<float>& sum, const std::vector<float>& sumsq,
                                 int length, int wlen) {
  std::vector<float> tstat(length + 1, 0.0f);
  const float eta = std::numeric_limits<float>::min();

  if (length < 2 * wlen || wlen < 2) return tstat;

  for (int i = wlen; i <= length - wlen; ++i) {
    float sum1 = sum[i];
    float sumsq1 = sumsq[i];
    if (i > wlen) {
      sum1 -= sum[i - wlen];
      sumsq1 -= sumsq[i - wlen];
    }
    float sum2 = sum[i + wlen] - sum[i];
    float sumsq2 = sumsq[i + wlen] - sumsq[i];
    float mean1 = sum1 / wlen;
    float mean2 = sum2 / wlen;
    float combined_var = (sumsq1 / wlen - mean1 * mean1 + sumsq2 / wlen - mean2 * mean2) / wlen;
    combined_var = std::max(combined_var, eta);
    const float delta_mean = mean2 - mean1;
    tstat[i] = std::abs(delta_mean) / std::sqrt(combined_var);
  }

  return tstat;
}

std::vector<int> short_long_peak_detector(const std::vector<float>& t1,
                                          const std::vector<float>& t2, int s_len,
                                          int window_length1, int window_length2, float threshold1,
                                          float threshold2, float peak_height) {
  std::vector<int> peaks;

  struct Detector {
    int DEF_PEAK_POS = -1;
    float DEF_PEAK_VAL = std::numeric_limits<float>::max();
    const std::vector<float>* sig;
    int s_len;
    float threshold;
    int window_length;
    int masked_to = 0;
    int peak_pos = -1;
    float peak_value = std::numeric_limits<float>::max();
    bool valid_peak = false;
  };

  Detector short_detector{
      -1, std::numeric_limits<float>::max(), &t1,  s_len, threshold1, window_length1, 0,
      -1, std::numeric_limits<float>::max(), false};
  Detector long_detector{
      -1, std::numeric_limits<float>::max(), &t2,  s_len, threshold2, window_length2, 0,
      -1, std::numeric_limits<float>::max(), false};

  Detector* detectors[2] = {&short_detector, &long_detector};

  for (int i = 0; i < s_len; ++i) {
    for (int k = 0; k < 2; ++k) {
      Detector* detector = detectors[k];
      if (detector->masked_to >= i) continue;

      float current_value = (*detector->sig)[i];

      if (detector->peak_pos == detector->DEF_PEAK_POS) {
        if (current_value < detector->peak_value) {
          detector->peak_value = current_value;
        } else if (current_value - detector->peak_value > peak_height) {
          detector->peak_value = current_value;
          detector->peak_pos = i;
        }
      } else {
        if (current_value > detector->peak_value) {
          detector->peak_value = current_value;
          detector->peak_pos = i;
        }
        // Dominate other detectors if we're going to fire
        if (detector->peak_value > detector->threshold) {
          for (int n_d = k + 1; n_d < 2; ++n_d) {
            detectors[n_d]->masked_to = detector->peak_pos + detectors[0]->window_length;
            detectors[n_d]->peak_pos = detectors[n_d]->DEF_PEAK_POS;
            detectors[n_d]->peak_value = detectors[n_d]->DEF_PEAK_VAL;
            detectors[n_d]->valid_peak = false;
          }
        }
        if (detector->peak_value - current_value > peak_height &&
            detector->peak_value > detector->threshold) {
          detector->valid_peak = true;
        }
        if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
          peaks.push_back(detector->peak_pos);
          detector->peak_pos = detector->DEF_PEAK_POS;
          detector->peak_value = current_value;
          detector->valid_peak = false;
        }
      }
    }
  }

  return peaks;
}

// IQR-filtered mean: filter outliers using interquartile range, then compute mean.
float iqr_filtered_mean(std::vector<float> segment) {
  if (segment.empty()) return 0.0f;

  std::sort(segment.begin(), segment.end());
  const size_t n = segment.size();
  const float q1 = segment[n / 4];
  const float q3 = segment[3 * n / 4];
  const float iqr = q3 - q1;
  const float lower_bound = q1 - iqr;
  const float upper_bound = q3 + iqr;

  float sum = 0.0f;
  int count = 0;
  for (const auto v : segment) {
    if (v >= lower_bound && v <= upper_bound) {
      sum += v;
      ++count;
    }
  }

  return count > 0 ? sum / count : 0.0f;
}

}  // namespace

StandardEventPipeline::StandardEventPipeline(EventPipelineConfig config)
    : config_(std::move(config)) {}

NormalizedSignal StandardEventPipeline::process(const io::RawRead& read) const {
  NormalizedSignal result;
  result.sampling_rate_hz = read.sampling_rate_hz;

  if (read.raw_signal.empty()) {
    return result;
  }

  // Step 1: Convert raw ADC values to picoamps
  const float raw_unit = (read.digitisation == 0.0f) ? 1.0f : (read.range / read.digitisation);
  std::vector<float> raw;
  raw.reserve(read.raw_signal.size());
  for (auto value : read.raw_signal) {
    raw.push_back((static_cast<float>(value) + read.offset) * raw_unit);
  }

  // Step 2: Normalize raw signal first
  // Compute running mean and stddev, then zscore normalize and filter outliers
  double sum = 0.0;
  double sum2 = 0.0;
  for (const auto v : raw) {
    sum += v;
    sum2 += v * v;
  }
  const double mean = sum / raw.size();
  const double stddev = std::sqrt(sum2 / raw.size() - mean * mean);

  std::vector<float> normalized;
  normalized.reserve(raw.size());
  for (const auto v : raw) {
    float norm_val = static_cast<float>((v - mean) / stddev);
    // Filter out samples where |z| >= 3
    if (norm_val > -3.0f && norm_val < 3.0f) {
      normalized.push_back(norm_val);
    }
  }

  if (normalized.empty()) {
    return result;
  }

  const int n_signals = static_cast<int>(normalized.size());

  // Step 3: Compute prefix sums on normalized signal
  std::vector<float> prefix_sum, prefix_sumsq;
  compute_sum_sumsq(normalized, prefix_sum, prefix_sumsq);

  // Step 4: Compute t-statistics
  auto tstat1 = compute_tstat(prefix_sum, prefix_sumsq, n_signals, config_.window_length1);
  auto tstat2 = compute_tstat(prefix_sum, prefix_sumsq, n_signals, config_.window_length2);

  // Step 5: Detect peaks
  auto peaks = short_long_peak_detector(tstat1, tstat2, n_signals, config_.window_length1,
                                        config_.window_length2, config_.threshold1,
                                        config_.threshold2, config_.peak_height);

  if (peaks.empty()) {
    return result;
  }

  // Step 6: Generate events using IQR-filtered mean (skip segments > 500 samples)
  constexpr int kMaxSegmentLength = 500;

  result.samples.reserve(peaks.size() + 1);

  int start_idx = 0;
  for (size_t pi = 0; pi < peaks.size(); ++pi) {
    if (peaks[pi] <= 0 || peaks[pi] >= n_signals) continue;

    int segment_length = peaks[pi] - start_idx;
    if (segment_length > 0 && segment_length < kMaxSegmentLength) {
      std::vector<float> segment(normalized.begin() + start_idx, normalized.begin() + peaks[pi]);
      result.samples.push_back(iqr_filtered_mean(segment));
    }
    start_idx = peaks[pi];
  }

  // Last segment
  int last_length = n_signals - start_idx;
  if (last_length > 0 && last_length < kMaxSegmentLength) {
    std::vector<float> segment(normalized.begin() + start_idx, normalized.end());
    result.samples.push_back(iqr_filtered_mean(segment));
  }

  // Note: Events are already normalized (from step 2), so no further normalization needed.
  // The IQR-filtered means are in z-score space.

  return result;
}

NormalizedSignal StandardEventPipeline::process_chunk(const float* pA, std::size_t len,
                                                     NormState& norm_state) const {
  NormalizedSignal result;
  if (len == 0) return result;

  /* Accumulate normalization stats (progressive) */
  for (std::size_t i = 0; i < len; ++i) {
    norm_state.sum += pA[i];
    norm_state.sum_sq += static_cast<double>(pA[i]) * pA[i];
  }
  norm_state.n += len;

  const double mean = norm_state.sum / static_cast<double>(norm_state.n);
  const double stddev =
      std::sqrt(norm_state.sum_sq / static_cast<double>(norm_state.n) - mean * mean);
  if (stddev < 1e-12) return result;

  /* Normalize this chunk and filter outliers */
  std::vector<float> normalized;
  normalized.reserve(len);
  for (std::size_t i = 0; i < len; ++i) {
    float norm_val = static_cast<float>((pA[i] - mean) / stddev);
    if (norm_val > -3.0f && norm_val < 3.0f) {
      normalized.push_back(norm_val);
    }
  }
  if (normalized.empty()) return result;

  // Stash normalized signal for debug dumping from caller
  result.debug_normalized = normalized;

  const int n_signals = static_cast<int>(normalized.size());

  /* Event detection (independent per chunk) */
  std::vector<float> prefix_sum, prefix_sumsq;
  compute_sum_sumsq(normalized, prefix_sum, prefix_sumsq);

  auto tstat1 = compute_tstat(prefix_sum, prefix_sumsq, n_signals, config_.window_length1);
  auto tstat2 = compute_tstat(prefix_sum, prefix_sumsq, n_signals, config_.window_length2);

  auto peaks = short_long_peak_detector(tstat1, tstat2, n_signals, config_.window_length1,
                                        config_.window_length2, config_.threshold1,
                                        config_.threshold2, config_.peak_height);
  if (peaks.empty()) return result;

  /* Generate events using IQR-filtered mean */
  constexpr int kMaxSegmentLength = 500;
  result.samples.reserve(peaks.size() + 1);

  int start_idx = 0;
  for (std::size_t pi = 0; pi < peaks.size(); ++pi) {
    if (peaks[pi] <= 0 || peaks[pi] >= n_signals) continue;
    int segment_length = peaks[pi] - start_idx;
    if (segment_length > 0 && segment_length < kMaxSegmentLength) {
      std::vector<float> segment(normalized.begin() + start_idx, normalized.begin() + peaks[pi]);
      result.samples.push_back(iqr_filtered_mean(segment));
    }
    start_idx = peaks[pi];
  }

  int last_length = n_signals - start_idx;
  if (last_length > 0 && last_length < kMaxSegmentLength) {
    std::vector<float> segment(normalized.begin() + start_idx, normalized.end());
    result.samples.push_back(iqr_filtered_mean(segment));
  }

  return result;
}

const EventPipelineConfig& StandardEventPipeline::config() const { return config_; }

std::string StandardEventPipeline::name() const { return "standard"; }

}  // namespace piru::signal
