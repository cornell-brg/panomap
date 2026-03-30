// SPDX-License-Identifier: MIT
// Interface for k-mer pore models.

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace piru::io {

class KmerModel {
public:
  virtual ~KmerModel() = default;

  // Human-readable model name (e.g., r10.4).
  virtual std::string name() const = 0;

  // k-mer length.
  virtual int k() const = 0;

  // Look up mean/current for a k-mer. Returns mean if found.
  virtual bool lookup(const std::string& kmer, double& mean) const = 0;

  // Build flat lookup table: float[4^k] indexed by 2-bit encoded k-mer integer.
  // Encoding: A=0, C=1, G=2, T=3. K-mer "ACG" (k=3) -> 0b00_01_10 = 6.
  // Entries for k-mers not in the model are set to 0.
  std::vector<float> buildFlatLookup() const;
};

}  // namespace piru::io
