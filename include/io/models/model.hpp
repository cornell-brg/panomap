// SPDX-License-Identifier: MIT
// Interface for k-mer pore models.

#pragma once

#include <string>
#include <utility>

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
};

}  // namespace piru::io
