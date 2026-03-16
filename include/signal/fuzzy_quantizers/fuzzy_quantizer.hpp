// SPDX-License-Identifier: MIT
// Interface for fuzzy quantization used during seeding.

#pragma once

#include <memory>
#include <string>

#include "signal/signal_types.hpp"

namespace piru::signal {

struct FuzzyQuantizerConfig {
  std::string backend{"piru"};  // Original default: "rawhash2"
  std::string params;
  std::string pore_model;  // For chemistry-specific defaults (r9/r10)
  float fine_min{-2.0f};
  float fine_max{2.0f};
  float fine_range{0.85f};  // Original default: 0.4f (DEV025: R9=0.9, R10=0.8)
  float diff{0.35f};         // Skip events within diff of last emitted (0 = disabled, RH2: 0.35)
  std::uint32_t qbits{4};
  std::uint32_t n_bins{10};  // Original default: 0 (=16). DEV025: 10 works well
  bool operator==(const FuzzyQuantizerConfig& other) const = default;
};

class FuzzyQuantizer {
public:
  virtual ~FuzzyQuantizer() = default;

  virtual FuzzyQuantizedSignal quantize(const NormalizedSignal& signal) const = 0;
  virtual const FuzzyQuantizerConfig& config() const = 0;
  virtual std::string name() const = 0;
};

using FuzzyQuantizerPtr = std::unique_ptr<FuzzyQuantizer>;

}  // namespace piru::signal
