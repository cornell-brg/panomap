// SPDX-License-Identifier: MIT

#include "index/path_walk_indexer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <limits>
#include <mutex>
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

PathWalkIndexResult pathWalkIndex(const AlnGraph& graph, const io::KmerModel& model,
                                  const signal::FuzzyQuantizer& fuzzy_quantizer,
                                  const signal::SeedExtractor& extractor,
                                  const PathWalkIndexConfig& config) {
  PathWalkIndexResult result;
  result.seed_store = std::make_unique<HashSeedStore>();
  result.linearization_coords.resize(graph.nodeCount());
  result.path_lengths.resize(graph.pathCount(), 0);

  const int pore_k = model.k();

  // Store extractor metadata for serialization
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

  // Optional: dump per-path normalization stats to file
  // NOTE: Disabled in parallel mode (would need mutex protection)
  std::ofstream norm_stats_file;
  if (!config.dump_norm_stats_path.empty()) {
    if (config.executor) {
      LOG_WARN("--dump-norm-stats disabled in parallel mode");
    } else {
      norm_stats_file.open(config.dump_norm_stats_path);
      if (norm_stats_file) {
        norm_stats_file << "path_name\tmean\tstddev\tnum_kmers\n";
      } else {
        LOG_WARN("Could not open norm stats file: " + config.dump_norm_stats_path);
      }
    }
  }

  // Atomic counters for parallel execution
  std::atomic<std::size_t> atomic_seeds_extracted{0};
  std::atomic<std::size_t> atomic_total_path_length{0};

  // Per-node mutexes for linearization_coords (multiple paths touch same nodes)
  std::vector<std::mutex> node_mutexes(graph.nodeCount());

  // Lambda to process a single path
  auto process_path = [&](std::size_t path_idx, HashSeedStore& thread_store) {
    const auto& path = graph.paths()[path_idx];

    // Step 1: Concatenate node sequences, track boundaries
    std::string path_sequence;
    std::vector<NodeBoundary> boundaries;

    for (const auto& step : path.steps) {
      std::size_t node_id = std::stoull(step.node_id);
      if (node_id >= graph.nodeCount()) continue;

      const auto& node = graph.node(node_id);
      boundaries.push_back({node_id, path_sequence.size()});
      path_sequence += node.sequence;
    }

    if (path_sequence.size() < static_cast<std::size_t>(pore_k)) {
      return;  // Path too short
    }

    // Step 2: Squigglize with single-pass variance (sum and sum_sq)
    std::vector<float> raw_signal;
    raw_signal.reserve(path_sequence.size() - pore_k + 1);

    double sum = 0.0;
    double sum_sq = 0.0;
    std::size_t count = 0;
    std::string kmer_buf(static_cast<std::size_t>(pore_k), '\0');

    for (std::size_t i = 0; i + static_cast<std::size_t>(pore_k) <= path_sequence.size(); ++i) {
      std::copy_n(path_sequence.data() + i, static_cast<std::size_t>(pore_k), kmer_buf.begin());

      if (hasNBase(kmer_buf)) {
        raw_signal.push_back(std::numeric_limits<float>::quiet_NaN());
        continue;
      }

      double mean = 0.0;
      if (!model.lookup(kmer_buf, mean)) {
        mean = 0.0;
      }
      float val = static_cast<float>(mean);
      raw_signal.push_back(val);
      sum += val;
      sum_sq += static_cast<double>(val) * static_cast<double>(val);
      ++count;
    }

    if (count == 0) return;

    // Step 3: Compute per-path normalization stats (single-pass)
    const double path_mean = sum / static_cast<double>(count);
    const double mean_sq = sum_sq / static_cast<double>(count);
    const double variance = mean_sq - (path_mean * path_mean);
    const double path_std = (variance > 0.0) ? std::sqrt(variance) : 0.0;

    // Step 4: Normalize and fuzzy quantize
    signal::NormalizedSignal normalized;
    normalized.samples.reserve(raw_signal.size());
    for (float val : raw_signal) {
      if (std::isnan(val)) {
        normalized.samples.push_back(val);
      } else {
        normalized.samples.push_back(safeNormalize(val, path_mean, path_std));
      }
    }

    auto fuzzy = fuzzy_quantizer.quantize(normalized);

    // Update path length (each path writes to its own slot)
    result.path_lengths[path_idx] = path_sequence.size();
    atomic_total_path_length.fetch_add(fuzzy.tokens.size(), std::memory_order_relaxed);

    // Step 5: Record linear coordinates (protected by per-node mutex)
    for (const auto& boundary : boundaries) {
      std::lock_guard<std::mutex> lock(node_mutexes[boundary.node_id]);
      result.linearization_coords[boundary.node_id].emplace_back(
          path_idx, static_cast<std::int64_t>(boundary.base_start));
    }

    // Step 6: Extract seeds into per-path store
    const auto seeds = extractor.extract(fuzzy);
    atomic_seeds_extracted.fetch_add(seeds.seeds.size(), std::memory_order_relaxed);

    HashSeedStore path_seed_store;
    for (const auto& seed : seeds.seeds) {
      auto [node_id, local_offset] = signalPositionToNodeOffset(seed.position, boundaries);
      path_seed_store.insert(seed.hash, SeedEntry{node_id, local_offset, seed.length});
    }

    // Step 7: Per-path frequency subsampling
    // - Below threshold (seed_freq_cutoff percentile): pass through
    // - Above threshold: subsample down to seed_freq_cap percentile freq
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
        subsample_cap = 0;  // harddrop mode
      }
    }

    // Merge seeds into thread-local store
    std::mt19937 rng(static_cast<unsigned>(path_idx));
    for (auto& [hash, hits] : path_seed_store.mutableData()) {
      if (hits.size() <= threshold) {
        // Below threshold -- merge all
        for (const auto& hit : hits) {
          thread_store.insert(hash, hit);
        }
      } else if (subsample_cap > 0) {
        // Above threshold -- subsample to subsample_cap
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
      // else: subsample_cap == 0 -> harddrop (discard entirely)
    }
  };

  const std::size_t num_paths = graph.pathCount();

#ifdef PIRU_USE_TBB
  if (config.executor) {
    // Parallel execution with thread-local stores
    tbb::enumerable_thread_specific<HashSeedStore> thread_stores;

    config.executor->parallel_for(0, num_paths, 1, [&](std::size_t path_idx) {
      process_path(path_idx, thread_stores.local());
    });

    // Merge thread-local stores into result
    for (auto& store : thread_stores) {
      result.seed_store->merge(store);
    }
  } else
#endif
  {
    // Sequential execution
    for (std::size_t path_idx = 0; path_idx < num_paths; ++path_idx) {
      process_path(path_idx, *result.seed_store);

      // Log per-path stats in sequential mode
      if (!config.executor) {
        const auto& path = graph.paths()[path_idx];
        LOG_DEBUG("Processed path " + path.name);
      }
    }
  }

  result.seeds_extracted = atomic_seeds_extracted.load();
  result.total_path_length = atomic_total_path_length.load();

  // Dedup seeds from shared regions across paths
  result.seed_store->deduplicate();
  result.seeds_unique = result.seed_store->size();

  // Compute max frequency
  std::size_t max_freq = 0;
  for (const auto& [hash, hits] : result.seed_store->data()) {
    max_freq = std::max(max_freq, hits.size());
  }
  result.seed_store->set_max_hash_frequency(max_freq);

  // Set frequency threshold to max_freq + 1 (no query-time filtering).
  // Per-path filtering already removed high-frequency seeds during indexing,
  // so we don't need a second round of filtering at query time.
  std::size_t threshold = max_freq + 1;
  result.seed_store->set_frequency_threshold(threshold);

  LOG_INFO("PathWalkIndex: " + std::to_string(result.seeds_extracted) + " seeds extracted, " +
           std::to_string(result.seeds_unique) + " unique (max_freq=" + std::to_string(max_freq) +
           ", global_threshold=" + std::to_string(threshold) + ")");

  return result;
}

}  // namespace piru::index
