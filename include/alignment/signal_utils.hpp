// SPDX-License-Identifier: MIT
// Signal utilities for alignment.

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <variant>

#include "signal/signal_types.hpp"

namespace piru::alignment {

/// Get the length of an AlignmentQuantizedSignal.
inline size_t signalLength(const signal::AlignmentQuantizedSignal& sig) {
  return std::visit([](const auto& vec) { return vec.size(); }, sig.data);
}

/// Slice an AlignmentQuantizedSignal from start to end (exclusive).
/// Returns a new signal with the same quantization kind.
inline signal::AlignmentQuantizedSignal sliceSignal(
    const signal::AlignmentQuantizedSignal& sig, size_t start, size_t end) {
  signal::AlignmentQuantizedSignal result;
  result.kind = sig.kind;

  std::visit(
      [&result, start, end](const auto& vec) {
        using T = typename std::decay_t<decltype(vec)>::value_type;
        size_t s = std::min(start, vec.size());
        size_t e = std::min(end, vec.size());
        if (s >= e) {
          result.data = std::vector<T>{};
          return;
        }
        result.data = std::vector<T>(vec.begin() + s, vec.begin() + e);
      },
      sig.data);

  return result;
}

/// Concatenate two AlignmentQuantizedSignals.
/// Both must have the same quantization kind.
inline signal::AlignmentQuantizedSignal concatSignals(
    const signal::AlignmentQuantizedSignal& a,
    const signal::AlignmentQuantizedSignal& b) {
  if (a.kind != b.kind) {
    throw std::invalid_argument("Cannot concatenate signals with different quantization kinds");
  }

  signal::AlignmentQuantizedSignal result;
  result.kind = a.kind;

  std::visit(
      [&result, &b](const auto& vec_a) {
        using T = typename std::decay_t<decltype(vec_a)>::value_type;
        const auto& vec_b = std::get<std::vector<T>>(b.data);
        std::vector<T> combined;
        combined.reserve(vec_a.size() + vec_b.size());
        combined.insert(combined.end(), vec_a.begin(), vec_a.end());
        combined.insert(combined.end(), vec_b.begin(), vec_b.end());
        result.data = std::move(combined);
      },
      a.data);

  return result;
}

/// Get signal value at index as float (for distance computation).
/// Works regardless of underlying quantization type.
inline float signalValueAt(const signal::AlignmentQuantizedSignal& sig, size_t idx) {
  return std::visit(
      [idx](const auto& vec) -> float {
        if (idx >= vec.size()) {
          throw std::out_of_range("Signal index out of range");
        }
        return static_cast<float>(vec[idx]);
      },
      sig.data);
}

/// Check if two signals have the same quantization kind.
inline bool sameQuantizationKind(const signal::AlignmentQuantizedSignal& a,
                                  const signal::AlignmentQuantizedSignal& b) {
  return a.kind == b.kind;
}

}  // namespace piru::alignment
