// SPDX-License-Identifier: MIT
// Landmark tokenizer implementation.
//
// Matches python prototype: pyrh2/peak_fingerprint.py
//   find_valleys() -> extract_peaks() -> _amp_peak_token()

#include "signal/tokenizers/landmark_tokenizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace panomap::signal {

namespace {

struct Peak {
  std::size_t left_idx;
  std::size_t apex_idx;
  std::size_t right_idx;
};

// Find local minima with prominence filtering.
// Matches python find_valleys() in peak_fingerprint.py:45-96
std::vector<std::size_t> find_valleys(const float* samples, std::size_t n,
                                      float min_prominence) {
  if (n < 3) return {};

  // Pass 1: collect all local minima
  std::vector<std::size_t> raw;
  for (std::size_t i = 1; i + 1 < n; ++i) {
    if (samples[i] <= samples[i - 1] && samples[i] <= samples[i + 1]) {
      raw.push_back(i);
    }
  }

  // Include boundary points if they are local lows
  if (n >= 2 && samples[0] <= samples[1]) {
    raw.insert(raw.begin(), 0);
  }
  if (n >= 2 && samples[n - 1] <= samples[n - 2]) {
    raw.push_back(n - 1);
  }

  if (raw.size() < 2) return raw;

  // Pass 2: merge valleys with too-small peaks between them
  std::vector<std::size_t> filtered;
  filtered.push_back(raw[0]);

  for (std::size_t i = 1; i < raw.size(); ++i) {
    std::size_t left = filtered.back();
    std::size_t right = raw[i];

    // Find apex between left and right
    float apex_val = samples[left];
    for (std::size_t j = left + 1; j <= right; ++j) {
      apex_val = std::max(apex_val, samples[j]);
    }

    float avg_valley = (samples[left] + samples[right]) / 2.0f;
    float prominence = apex_val - avg_valley;

    if (prominence >= min_prominence) {
      filtered.push_back(right);
    } else {
      // Small peak: keep the deeper valley
      if (samples[right] < samples[left]) {
        filtered.back() = right;
      }
    }
  }

  return filtered;
}

// Extract peaks between consecutive valleys.
// Matches python extract_peaks() in peak_fingerprint.py:99-128
std::vector<Peak> extract_peaks(const float* samples, const std::vector<std::size_t>& valleys) {
  std::vector<Peak> peaks;
  if (valleys.size() < 2) return peaks;

  peaks.reserve(valleys.size() - 1);
  for (std::size_t i = 0; i + 1 < valleys.size(); ++i) {
    std::size_t left = valleys[i];
    std::size_t right = valleys[i + 1];

    // Find apex in [left, right]
    std::size_t apex = left;
    float apex_val = samples[left];
    for (std::size_t j = left + 1; j <= right; ++j) {
      if (samples[j] > apex_val) {
        apex_val = samples[j];
        apex = j;
      }
    }

    peaks.push_back({left, apex, right});
  }

  return peaks;
}

// Log-scale quantization: finer near 0, coarser for large values.
// Matches python _log_quantize() in peak_fingerprint.py:1092-1097
std::uint32_t log_quantize(float val, std::uint32_t bits, float max_val) {
  std::uint32_t n = 1u << bits;
  val = std::clamp(std::abs(val), 0.0f, max_val);
  float norm = std::log1p(val) / std::log1p(max_val);
  return std::min(static_cast<std::uint32_t>(norm * n), n - 1);
}

}  // namespace

LandmarkTokenizer::LandmarkTokenizer(TokenizerConfig config, LandmarkTokenizerConfig landmark_config)
    : config_(std::move(config)), landmark_config_(landmark_config) {}

TokenizedSignal LandmarkTokenizer::quantize(const NormalizedSignal& signal) const {
  TokenizedSignal result;
  const auto& samples = signal.samples;
  const std::size_t n = samples.size();
  if (n < 3) return result;

  auto valleys = find_valleys(samples.data(), n, landmark_config_.min_prominence);
  auto peaks = extract_peaks(samples.data(), valleys);

  if (peaks.empty()) return result;

  result.tokens.reserve(peaks.size());
  result.original_positions.reserve(peaks.size());

  const auto& input_pos = signal.original_positions;
  const bool has_input_pos = !input_pos.empty();

  for (const auto& peak : peaks) {
    float rise = samples[peak.apex_idx] - samples[peak.left_idx];
    float fall = samples[peak.apex_idx] - samples[peak.right_idx];

    auto rb = log_quantize(rise, landmark_config_.rise_bits, landmark_config_.max_amp);
    auto fb = log_quantize(fall, landmark_config_.fall_bits, landmark_config_.max_amp);
    auto token = static_cast<std::int16_t>((rb << landmark_config_.fall_bits) | fb);

    result.tokens.push_back(token);
    // Map through diff filter positions if present, otherwise use directly
    auto idx = static_cast<std::uint32_t>(peak.right_idx);
    result.original_positions.push_back(has_input_pos ? input_pos[idx] : idx);
  }

  return result;
}

const TokenizerConfig& LandmarkTokenizer::config() const { return config_; }

std::string LandmarkTokenizer::name() const { return "landmark"; }

}  // namespace panomap::signal
