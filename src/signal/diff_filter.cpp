// SPDX-License-Identifier: MIT

#include "signal/diff_filter.hpp"

#include <cmath>
#include <limits>

namespace panomap::signal {

void apply_diff_filter(NormalizedSignal& signal, float diff) {
  if (diff <= 0.0f) return;

  std::vector<float> filtered;
  std::vector<std::uint32_t> positions;
  filtered.reserve(signal.samples.size());
  positions.reserve(signal.samples.size());

  float last_emitted = std::numeric_limits<float>::quiet_NaN();

  for (std::uint32_t i = 0; i < signal.samples.size(); ++i) {
    const float sample = signal.samples[i];

    // Always keep NaN sentinels
    if (std::isnan(sample)) {
      filtered.push_back(sample);
      positions.push_back(i);
      continue;
    }

    // Drop events too similar to last emitted
    if (!std::isnan(last_emitted) && std::fabs(sample - last_emitted) < diff) {
      continue;
    }

    last_emitted = sample;
    filtered.push_back(sample);
    positions.push_back(i);
  }

  signal.samples = std::move(filtered);
  signal.original_positions = std::move(positions);
}

}  // namespace panomap::signal
