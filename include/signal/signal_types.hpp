// SPDX-License-Identifier: MIT
// Shared signal-processing data structures used across mapping and indexing.

#pragma once

#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

namespace piru::signal {

struct NormalizedSignal {
    std::vector<float> samples;
    float sampling_rate_hz{0.0f};
};

struct SignalEvent {
    std::size_t start{0};
    std::size_t length{0};
    float mean{0.0f};
    float stdv{0.0f};
};

struct EventSeries {
    std::vector<SignalEvent> events;
};

struct FuzzyQuantizedSignal {
    std::vector<std::int16_t> tokens;
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
    std::size_t length{0};     // Seed coverage length (initially k, may increase after merging)
};

struct SeedBuffer {
    std::vector<Seed> seeds;
};

}  // namespace piru::signal
