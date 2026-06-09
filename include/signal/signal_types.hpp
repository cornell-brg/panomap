// SPDX-License-Identifier: MIT
// Shared signal-processing data structures used across mapping and indexing.

#pragma once

#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

namespace panomap::signal {

struct NormalizedSignal {
  std::vector<float> samples;
  float sampling_rate_hz{0.0f};
  std::vector<float> debug_normalized;  // DEV-62: pre-event-detection normalized signal for tracing
  // When diff filter compresses the event stream, this maps compressed sample
  // index -> original event index. Empty when no compression (positions = indices).
  std::vector<std::uint32_t> original_positions;
};

struct SignalEvent {
  std::size_t start{0};
  std::size_t length{0};
  float mean{0.0f};
  float stdv{0.0f};
};

struct EventSeries {
  std::vector<SignalEvent> events;
  float sampling_rate_hz{0.0f};
};

struct TokenizedSignal {
  std::vector<std::int16_t> tokens;
  // Passed through from NormalizedSignal. Maps compressed token index -> original
  // event index. Empty when no compression (positions = indices).
  std::vector<std::uint32_t> original_positions;
};

enum class AlignmentQuantizationKind { kFloat32, kInt16, kInt8 };

using AlignmentQuantizedPayload =
    std::variant<std::vector<float>, std::vector<std::int16_t>, std::vector<std::int8_t>>;

struct AlignmentQuantizedSignal {
  AlignmentQuantizationKind kind{AlignmentQuantizationKind::kFloat32};
  AlignmentQuantizedPayload data;
};

struct Seed {
  std::uint64_t hash{0};
  std::size_t position{0};
  std::size_t length{0};  // Seed coverage length (initially k, may increase after merging)
};

struct SeedBuffer {
  std::vector<Seed> seeds;
};

}  // namespace panomap::signal
