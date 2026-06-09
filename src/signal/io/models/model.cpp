// SPDX-License-Identifier: MIT

#include "signal/io/models/model.hpp"

namespace panomap::io {

std::vector<float> KmerModel::buildFlatLookup() const {
  const int kk = k();
  const std::size_t table_size = 1ULL << (2 * kk);  // 4^k
  std::vector<float> flat(table_size, 0.0f);

  // Enumerate all 4^k k-mers, look up each in the hash table
  std::string kmer(static_cast<std::size_t>(kk), 'A');
  constexpr char bases[] = {'A', 'C', 'G', 'T'};

  for (std::size_t idx = 0; idx < table_size; ++idx) {
    // Convert integer index to k-mer string
    std::size_t val = idx;
    for (int j = kk - 1; j >= 0; --j) {
      kmer[static_cast<std::size_t>(j)] = bases[val & 3];
      val >>= 2;
    }

    double mean = 0.0;
    if (lookup(kmer, mean)) {
      flat[idx] = static_cast<float>(mean);
    }
  }

  return flat;
}

}  // namespace panomap::io
