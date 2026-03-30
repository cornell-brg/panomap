// SPDX-License-Identifier: MIT

#include "index/path_walk_indexer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <limits>
// mutex no longer needed -- linearization uses thread-local merge
#include <random>

#ifdef PIRU_USE_TBB
#include <tbb/enumerable_thread_specific.h>
#endif

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::index {

namespace {

/* Check if k-mer contains N or other non-ACGT bases. */
bool hasNBase(const std::string_view& kmer) {
  for (char c : kmer) {
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      return true;
    }
  }
  return false;
}

/* Z-score normalization with outlier clipping to [-3, 3]. */
float safeNormalize(float value, double mean, double stddev) {
  if (stddev < 1e-6) return 0.0f;
  float z_score = static_cast<float>((value - mean) / stddev);
  return std::clamp(z_score, -3.0f, 3.0f);
}

/* Track where each node starts in the concatenated path sequence. */
struct NodeBoundary {
  std::size_t node_id;
  std::size_t base_start;
};

/* Map signal position back to (node_id, local_offset). Binary search since boundaries are sorted.
 */
std::pair<std::size_t, std::size_t> signalPositionToNodeOffset(
    std::size_t signal_pos, const std::vector<NodeBoundary>& boundaries) {
  // Find first boundary with base_start > signal_pos, then back up one
  auto it =
      std::upper_bound(boundaries.begin(), boundaries.end(), signal_pos,
                       [](std::size_t pos, const NodeBoundary& b) { return pos < b.base_start; });
  if (it != boundaries.begin()) {
    --it;
  }
  return {it->node_id, signal_pos - it->base_start};
}

}  // namespace

// --- FlatGraph implementation ---

PathWalkIndexResult pathWalkIndex(const FlatGraph& graph, const io::KmerModel& model,
                                  const signal::FuzzyQuantizer& fuzzy_quantizer,
                                  const signal::SeedExtractor& extractor,
                                  const PathWalkIndexConfig& config) {
  PathWalkIndexResult result;
  result.seed_store = std::make_unique<HashSeedStore>();
  result.linearization_coords.resize(graph.nodeCount());
  result.path_lengths.resize(graph.pathCount(), 0);

  const int pore_k = model.k();
  const auto pore_flat = model.buildFlatLookup();  // float[4^k], built once
  const std::uint64_t kmer_mask = (1ULL << (2 * pore_k)) - 1;

  result.seed_store->set_extractor_name(extractor.name());
  std::map<std::string, std::string> params;
  const auto& cfg = extractor.config();
  params["backend"] = cfg.backend;
  params["k"] = std::to_string(cfg.k);
  params["stride"] = std::to_string(cfg.stride);
  params["qbits"] = std::to_string(cfg.qbits);
  params["window"] = std::to_string(cfg.window);
  result.seed_store->set_params(std::move(params));
  result.seed_store->set_filter_fraction(config.seed_freq_cutoff);

  std::atomic<std::size_t> atomic_seeds_extracted{0};
  std::atomic<std::size_t> atomic_total_path_length{0};

  using LinearCoordVec = std::vector<std::vector<index::LinearCoordinate>>;

  auto process_path = [&](std::size_t path_idx, HashSeedStore& thread_store,
                          LinearCoordVec& thread_coords) {
    // Step 1: Compute path length and build boundaries
    std::vector<NodeBoundary> boundaries;

    const auto* steps = graph.pathStepsBegin(static_cast<std::uint32_t>(path_idx));
    std::size_t num_steps = graph.pathStepCount(static_cast<std::uint32_t>(path_idx));

    std::size_t path_len = 0;
    for (std::size_t s = 0; s < num_steps; ++s) {
      std::uint32_t node_id = steps[s];
      if (node_id >= graph.nodeCount()) continue;
      boundaries.push_back({node_id, path_len});
      path_len += graph.seqLen(node_id);
    }

    if (path_len < static_cast<std::size_t>(pore_k)) return;

    // Step 2: Fast squigglize using rolling 2-bit k-mer + flat pore model array.
    // Iterates through node sequences directly from 2-bit packed storage.
    // No string allocation, no hash table lookup.
    const std::size_t num_positions = path_len - static_cast<std::size_t>(pore_k) + 1;
    std::vector<float> raw_signal(num_positions);

    std::uint64_t kmer = 0;
    int valid_bases = 0;
    std::size_t sig_idx = 0;
    std::size_t base_pos = 0;

    for (std::size_t s = 0; s < num_steps; ++s) {
      std::uint32_t node_id = steps[s];
      if (node_id >= graph.nodeCount()) continue;
      std::uint32_t node_len = static_cast<std::uint32_t>(graph.seqLen(node_id));

      for (std::uint32_t j = 0; j < node_len; ++j) {
        if (graph.isN(node_id, j)) {
          // N base: reset rolling k-mer, emit NaN for affected positions
          kmer = 0;
          valid_bases = 0;
          if (base_pos >= static_cast<std::size_t>(pore_k) - 1) {
            // This N invalidates the current position
          }
          // Positions where the k-mer window includes this N will get NaN
          // We handle this by only writing valid values below
        } else {
          std::uint8_t base = graph.base2bit(node_id, j);
          kmer = ((kmer << 2) | base) & kmer_mask;
          valid_bases++;
        }

        if (base_pos >= static_cast<std::size_t>(pore_k) - 1) {
          if (valid_bases >= pore_k) {
            raw_signal[sig_idx] = pore_flat[kmer];
          } else {
            raw_signal[sig_idx] = std::numeric_limits<float>::quiet_NaN();
          }
          sig_idx++;
        }
        base_pos++;
      }
    }

    // Check if any valid signals were produced
    std::size_t count = 0;
    for (std::size_t i = 0; i < num_positions; ++i) {
      if (!std::isnan(raw_signal[i])) count++;
    }
    if (count == 0) return;

    signal::NormalizedSignal normalized;
    normalized.samples = std::move(raw_signal);

    auto fuzzy = fuzzy_quantizer.quantize(normalized);

    result.path_lengths[path_idx] = path_len;
    atomic_total_path_length.fetch_add(fuzzy.tokens.size(), std::memory_order_relaxed);

    // Step 5: Linear coordinates
    for (const auto& boundary : boundaries) {
      thread_coords[boundary.node_id].emplace_back(
          path_idx, static_cast<std::int64_t>(boundary.base_start));
    }

    // Step 6: Extract seeds
    const auto seeds = extractor.extract(fuzzy);
    atomic_seeds_extracted.fetch_add(seeds.seeds.size(), std::memory_order_relaxed);

    HashSeedStore path_seed_store;
    for (const auto& seed : seeds.seeds) {
      auto [node_id, local_offset] = signalPositionToNodeOffset(seed.position, boundaries);
      path_seed_store.insert(seed.hash, SeedEntry{
          static_cast<std::uint32_t>(node_id), static_cast<std::uint32_t>(local_offset)});
    }

    // Step 7: Frequency subsampling
    std::size_t threshold = std::numeric_limits<std::size_t>::max();
    std::size_t subsample_cap = std::numeric_limits<std::size_t>::max();
    if (!path_seed_store.data().empty() && config.seed_freq_cutoff < 1.0) {
      std::vector<std::size_t> path_frequencies;
      path_frequencies.reserve(path_seed_store.data().size());
      for (const auto& [hash, hits] : path_seed_store.data()) {
        path_frequencies.push_back(hits.size());
      }
      std::sort(path_frequencies.begin(), path_frequencies.end());

      auto freq_at = [&](double p) -> std::size_t {
        std::size_t pos = static_cast<std::size_t>(path_frequencies.size() * p);
        pos = std::min(pos, path_frequencies.size() - 1);
        return path_frequencies[pos];
      };
      threshold = freq_at(std::clamp(config.seed_freq_cutoff, 0.0, 1.0));
      if (config.seed_freq_cap > 0.0) {
        subsample_cap =
            std::max<std::size_t>(1, freq_at(std::clamp(config.seed_freq_cap, 0.0, 1.0)));
      } else {
        subsample_cap = 0;
      }
    }

    std::mt19937 rng(static_cast<unsigned>(path_idx));
    for (auto& [hash, hits] : path_seed_store.mutableData()) {
      if (hits.size() <= threshold) {
        for (const auto& hit : hits) {
          thread_store.insert(hash, hit);
        }
      } else if (subsample_cap > 0) {
        std::size_t target = std::min(subsample_cap, hits.size());
        std::size_t n = hits.size();
        for (std::size_t i = 0; i < target; ++i) {
          std::size_t j = i + (rng() % (n - i));
          std::swap(hits[i], hits[j]);
        }
        for (std::size_t i = 0; i < target; ++i) {
          thread_store.insert(hash, hits[i]);
        }
      }
    }
  };

  const std::size_t num_paths = graph.pathCount();

#ifdef PIRU_USE_TBB
  if (config.executor) {
    tbb::enumerable_thread_specific<HashSeedStore> thread_stores;
    tbb::enumerable_thread_specific<LinearCoordVec> thread_coords(
        [&]() { return LinearCoordVec(graph.nodeCount()); });

    config.executor->parallel_for(0, num_paths, 1, [&](std::size_t path_idx) {
      process_path(path_idx, thread_stores.local(), thread_coords.local());
    });

    for (auto& store : thread_stores) {
      result.seed_store->merge(store);
    }
    for (auto& coords : thread_coords) {
      for (std::size_t i = 0; i < coords.size(); ++i) {
        auto& dst = result.linearization_coords[i];
        dst.insert(dst.end(), coords[i].begin(), coords[i].end());
      }
    }
  } else
#endif
  {
    for (std::size_t path_idx = 0; path_idx < num_paths; ++path_idx) {
      process_path(path_idx, *result.seed_store, result.linearization_coords);
    }
  }

  result.seeds_extracted = atomic_seeds_extracted.load();
  result.total_path_length = atomic_total_path_length.load();

  result.seed_store->deduplicate();
  result.seeds_unique = result.seed_store->size();

  std::size_t max_freq = 0;
  for (const auto& [hash, hits] : result.seed_store->data()) {
    max_freq = std::max(max_freq, hits.size());
  }
  result.seed_store->set_max_hash_frequency(max_freq);

  std::size_t threshold = max_freq + 1;
  result.seed_store->set_frequency_threshold(threshold);

  LOG_INFO("PathWalkIndex [FlatGraph]: " + std::to_string(result.seeds_extracted) +
           " seeds extracted, " + std::to_string(result.seeds_unique) +
           " unique (max_freq=" + std::to_string(max_freq) + ")");

  return result;
}

}  // namespace piru::index
