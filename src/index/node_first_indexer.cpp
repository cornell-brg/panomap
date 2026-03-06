// SPDX-License-Identifier: MIT

#include "index/node_first_indexer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>

#ifdef PIRU_USE_TBB
#include <tbb/enumerable_thread_specific.h>
#endif

#include "signal/signal_types.hpp"
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
float safeNormalize(float value, float mean, float stddev) {
  if (stddev < 1e-6f) return 0.0f;
  float z_score = (value - mean) / stddev;
  return std::clamp(z_score, -3.0f, 3.0f);
}

/**
 * Get fuzzy samples in a hash window starting at a given position.
 *
 * This function computes fresh from path context:
 * 1. Collect sequence context starting at pos, extending into subsequent nodes
 * 2. Squigglize the needed k-mers
 * 3. Normalize with global mean/std
 * 4. Fuzzy quantize
 * 5. Return w fuzzy samples for hashing
 *
 * @param graph The alignment graph
 * @param model Pore model for squigglization
 * @param fuzzy_quantizer For converting normalized signal to tokens
 * @param path_idx Index of the path being walked
 * @param node_idx Index of the node within the path's steps
 * @param offset Position within the node to start the hash window
 * @param pore_k Pore model k-mer size
 * @param window_size Number of fuzzy samples needed (seed_k)
 * @param global_mean Global mean for normalization
 * @param global_std Global std for normalization
 *
 * @return Vector of fuzzy tokens (size = window_size if successful, empty if insufficient context)
 */
std::vector<std::int16_t> getHashWindow(const AlnGraph& graph, const io::KmerModel& model,
                                        const signal::FuzzyQuantizer& fuzzy_quantizer,
                                        std::size_t path_idx, std::size_t node_idx,
                                        std::size_t offset, int pore_k, std::size_t window_size,
                                        float global_mean, float global_std) {
  const auto& path = graph.paths()[path_idx];
  const std::size_t num_steps = path.steps.size();

  // We need (window_size) fuzzy samples
  // Each fuzzy sample requires (pore_k) bases
  // So we need (window_size + pore_k - 1) bases total starting from offset
  const std::size_t bases_needed = window_size + static_cast<std::size_t>(pore_k) - 1;

  // Collect bases from current node and subsequent nodes
  std::string context;
  context.reserve(bases_needed);

  std::size_t current_node_idx = node_idx;
  std::size_t current_offset = offset;

  while (context.size() < bases_needed && current_node_idx < num_steps) {
    std::size_t node_id = std::stoull(path.steps[current_node_idx].node_id);
    if (node_id >= graph.nodeCount()) {
      ++current_node_idx;
      current_offset = 0;
      continue;
    }

    const auto& node = graph.node(node_id);
    const std::size_t node_len = node.sequence.size();

    // How many bases can we take from this node?
    std::size_t available = (current_offset < node_len) ? (node_len - current_offset) : 0;
    std::size_t to_take = std::min(available, bases_needed - context.size());

    if (to_take > 0) {
      context.append(node.sequence, current_offset, to_take);
    }

    ++current_node_idx;
    current_offset = 0;
  }

  if (context.size() < bases_needed) {
    // Not enough context (end of path)
    return {};
  }

  // Squigglize: produce window_size k-mer values
  std::vector<float> raw_values;
  raw_values.reserve(window_size);
  std::string kmer_buf(static_cast<std::size_t>(pore_k), '\0');

  for (std::size_t i = 0; i < window_size; ++i) {
    std::copy_n(context.data() + i, static_cast<std::size_t>(pore_k), kmer_buf.begin());

    if (hasNBase(kmer_buf)) {
      // Contains N - can't hash this window
      return {};
    }

    double mean = 0.0;
    if (!model.lookup(kmer_buf, mean)) {
      mean = 0.0;
    }
    raw_values.push_back(static_cast<float>(mean));
  }

  // Normalize with global stats
  signal::NormalizedSignal normalized;
  normalized.samples.reserve(window_size);
  for (float val : raw_values) {
    normalized.samples.push_back(safeNormalize(val, global_mean, global_std));
  }

  // Fuzzy quantize
  auto fuzzy = fuzzy_quantizer.quantize(normalized);
  return std::move(fuzzy.tokens);
}

}  // namespace

NodeFirstIndexResult nodeFirstIndex(const AlnGraph& graph, const io::KmerModel& model,
                                    const signal::FuzzyQuantizer& fuzzy_quantizer,
                                    const signal::SeedExtractor& extractor,
                                    const NodeFirstIndexConfig& config) {
  NodeFirstIndexResult result;
  result.seed_store = std::make_unique<HashSeedStore>();
  result.linearization_coords.resize(graph.nodeCount());
  result.path_lengths.resize(graph.pathCount(), 0);

  const int pore_k = model.k();
  const std::size_t seed_k = config.seed_k;  // Hash window size

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

  // =========================================================================
  // Pass 1a: Squigglize node interiors, accumulate global stats
  // =========================================================================
  timing::start("nodefirst:squigglize");

  // Per-node raw k-mer values (temporary storage for pass 1)
  std::vector<std::vector<float>> node_raw_values(graph.nodeCount());

  // Single-pass stats: track sum and sum_sq for variance calculation
  // variance = (sum_sq / count) - (mean * mean)
  std::atomic<double> atomic_sum{0.0};
  std::atomic<double> atomic_sum_sq{0.0};
  std::atomic<std::size_t> atomic_count{0};

  const std::size_t num_nodes = graph.nodeCount();
  const std::size_t grain = 100;  // Process nodes in batches

  auto squigglize_node = [&](std::size_t node_id) {
    const auto& node = graph.node(node_id);
    const std::size_t node_len = node.sequence.size();

    // Can only squigglize positions [0, node_len - pore_k]
    if (node_len < static_cast<std::size_t>(pore_k)) {
      return;  // Node too short to produce any k-mer values
    }

    const std::size_t num_kmers = node_len - static_cast<std::size_t>(pore_k) + 1;
    auto& raw_values = node_raw_values[node_id];
    raw_values.reserve(num_kmers);

    std::string kmer_buf(static_cast<std::size_t>(pore_k), '\0');

    // Accumulate locally, then do single atomic update
    double local_sum = 0.0;
    double local_sum_sq = 0.0;
    std::size_t local_count = 0;

    for (std::size_t i = 0; i < num_kmers; ++i) {
      std::copy_n(node.sequence.data() + i, static_cast<std::size_t>(pore_k), kmer_buf.begin());

      if (hasNBase(kmer_buf)) {
        raw_values.push_back(std::numeric_limits<float>::quiet_NaN());
        continue;
      }

      double mean = 0.0;
      if (!model.lookup(kmer_buf, mean)) {
        mean = 0.0;
      }
      float val = static_cast<float>(mean);
      raw_values.push_back(val);
      local_sum += val;
      local_sum_sq += static_cast<double>(val) * static_cast<double>(val);
      ++local_count;
    }

    // Single atomic update per node (not per k-mer)
    if (local_count > 0) {
      atomic_sum.fetch_add(local_sum, std::memory_order_relaxed);
      atomic_sum_sq.fetch_add(local_sum_sq, std::memory_order_relaxed);
      atomic_count.fetch_add(local_count, std::memory_order_relaxed);
    }
  };

  if (config.executor) {
    // Parallel execution
    config.executor->parallel_for(0, num_nodes, grain, squigglize_node);
  } else {
    // Sequential fallback
    for (std::size_t node_id = 0; node_id < num_nodes; ++node_id) {
      squigglize_node(node_id);
    }
  }

  const double global_sum = atomic_sum.load();
  const double global_sum_sq = atomic_sum_sq.load();
  const std::size_t global_count = atomic_count.load();

  timing::stop("nodefirst:squigglize");

  // =========================================================================
  // Pass 1b: Compute global mean/std (single-pass from sum and sum_sq)
  // =========================================================================
  timing::start("nodefirst:norm_quant");

  if (global_count == 0) {
    LOG_WARN("NodeFirstIndex: No k-mer values found in graph");
    return result;
  }

  // Single-pass variance: Var(X) = E[X²] - E[X]²
  const double mean = global_sum / static_cast<double>(global_count);
  const double mean_sq = global_sum_sq / static_cast<double>(global_count);
  const double variance = mean_sq - (mean * mean);

  result.global_mean = static_cast<float>(mean);
  result.global_std = (variance > 0.0) ? static_cast<float>(std::sqrt(variance)) : 1.0f;

  LOG_DEBUG("NodeFirstIndex: global_mean=" + std::to_string(result.global_mean) +
            ", global_std=" + std::to_string(result.global_std) + " (from " +
            std::to_string(global_count) + " k-mers)");

  // =========================================================================
  // Pass 1c: Normalize, fuzzy quantize, index node interiors
  // =========================================================================

  std::atomic<std::size_t> atomic_seeds_interior{0};

  auto index_node_interior = [&](std::size_t node_id, HashSeedStore& local_store) {
    const auto& raw_values = node_raw_values[node_id];
    if (raw_values.empty()) return;

    // Normalize
    signal::NormalizedSignal normalized;
    normalized.samples.reserve(raw_values.size());
    for (float val : raw_values) {
      if (std::isnan(val)) {
        normalized.samples.push_back(val);
      } else {
        normalized.samples.push_back(safeNormalize(val, result.global_mean, result.global_std));
      }
    }

    // Fuzzy quantize
    auto fuzzy = fuzzy_quantizer.quantize(normalized);

    // Index positions with complete hash windows: [0, num_fuzzy - seed_k]
    // Position i can be hashed if we have fuzzy samples [i, i + seed_k)
    std::size_t local_count = 0;
    if (fuzzy.tokens.size() >= seed_k) {
      const std::size_t indexable_end = fuzzy.tokens.size() - seed_k + 1;

      for (std::size_t pos = 0; pos < indexable_end; pos += config.seed_stride) {
        signal::FuzzyQuantizedSignal window_fuzzy;
        window_fuzzy.tokens.assign(
            fuzzy.tokens.begin() + static_cast<std::ptrdiff_t>(pos),
            fuzzy.tokens.begin() + static_cast<std::ptrdiff_t>(pos + seed_k));

        auto seeds = extractor.extract(window_fuzzy);
        for (const auto& seed : seeds.seeds) {
          local_store.insert(seed.hash, SeedEntry{node_id, pos + seed.position, seed.length});
          ++local_count;
        }
      }
    }
    if (local_count > 0) {
      atomic_seeds_interior.fetch_add(local_count, std::memory_order_relaxed);
    }
  };

#ifdef PIRU_USE_TBB
  if (config.executor) {
    // Parallel execution with thread-local stores
    tbb::enumerable_thread_specific<HashSeedStore> interior_stores;

    config.executor->parallel_for(0, num_nodes, grain, [&](std::size_t node_id) {
      index_node_interior(node_id, interior_stores.local());
    });

    // Merge thread-local stores into result
    for (auto& store : interior_stores) {
      result.seed_store->merge(store);
    }
  } else
#endif
  {
    // Sequential fallback
    for (std::size_t node_id = 0; node_id < num_nodes; ++node_id) {
      index_node_interior(node_id, *result.seed_store);
    }
  }

  result.seeds_interior = atomic_seeds_interior.load();
  timing::stop("nodefirst:norm_quant");

  // Free raw values - no longer needed
  node_raw_values.clear();
  node_raw_values.shrink_to_fit();

  // =========================================================================
  // Pass 2: Boundary fill (path walks)
  // =========================================================================
  timing::start("nodefirst:hash_boundary");

  std::atomic<std::size_t> atomic_seeds_boundary{0};

  // Per-node mutexes for linearization_coords updates
  // (Multiple paths may touch the same node)
  std::vector<std::mutex> node_mutexes(graph.nodeCount());

  auto process_path = [&](std::size_t path_idx, HashSeedStore& local_store) {
    const auto& path = graph.paths()[path_idx];

    std::size_t path_base_pos = 0;  // Running base position in linearized path
    std::size_t local_seed_count = 0;

    for (std::size_t step_idx = 0; step_idx < path.steps.size(); ++step_idx) {
      std::size_t node_id = std::stoull(path.steps[step_idx].node_id);
      if (node_id >= graph.nodeCount()) continue;

      const auto& node = graph.node(node_id);
      const std::size_t node_len = node.sequence.size();

      // Record linear coordinate for this node on this path (base position)
      // Protected by per-node mutex since multiple paths may touch same node
      {
        std::lock_guard<std::mutex> lock(node_mutexes[node_id]);
        result.linearization_coords[node_id].emplace_back(path_idx,
                                                          static_cast<std::int64_t>(path_base_pos));
      }

      // Compute where interior indexing stopped
      std::size_t first_boundary_pos = 0;
      if (node_len >= static_cast<std::size_t>(pore_k) + seed_k - 1) {
        std::size_t num_fuzzy = node_len - static_cast<std::size_t>(pore_k) + 1;
        if (num_fuzzy >= seed_k) {
          first_boundary_pos = num_fuzzy - seed_k + 1;
        }
      }

      // Fill boundary positions using getHashWindow
      for (std::size_t pos = first_boundary_pos; pos < node_len; pos += config.seed_stride) {
        auto window = getHashWindow(graph, model, fuzzy_quantizer, path_idx, step_idx, pos, pore_k,
                                    seed_k, result.global_mean, result.global_std);

        if (window.empty()) {
          continue;
        }

        signal::FuzzyQuantizedSignal fuzzy_window;
        fuzzy_window.tokens = std::move(window);

        auto seeds = extractor.extract(fuzzy_window);
        for (const auto& seed : seeds.seeds) {
          local_store.insert(seed.hash, SeedEntry{node_id, pos + seed.position, seed.length});
          ++local_seed_count;
        }
      }

      path_base_pos += node_len;
    }

    // Record path length (each path writes to its own slot - no sync needed)
    result.path_lengths[path_idx] = path_base_pos;

    if (local_seed_count > 0) {
      atomic_seeds_boundary.fetch_add(local_seed_count, std::memory_order_relaxed);
    }
  };

  const std::size_t num_paths = graph.pathCount();

#ifdef PIRU_USE_TBB
  if (config.executor) {
    // Parallel execution with thread-local stores
    tbb::enumerable_thread_specific<HashSeedStore> boundary_stores;

    config.executor->parallel_for(0, num_paths, 1, [&](std::size_t path_idx) {
      process_path(path_idx, boundary_stores.local());
    });

    // Merge thread-local stores into result
    for (auto& store : boundary_stores) {
      result.seed_store->merge(store);
    }
  } else
#endif
  {
    // Sequential fallback
    for (std::size_t path_idx = 0; path_idx < num_paths; ++path_idx) {
      process_path(path_idx, *result.seed_store);
    }
  }

  result.seeds_boundary = atomic_seeds_boundary.load();
  timing::stop("nodefirst:hash_boundary");

  // =========================================================================
  // Finalize: dedup, compute stats
  // =========================================================================

  result.seed_store->deduplicate();
  result.seeds_unique = result.seed_store->size();

  // Compute max frequency
  std::size_t max_freq = 0;
  for (const auto& [hash, hits] : result.seed_store->data()) {
    max_freq = std::max(max_freq, hits.size());
  }
  result.seed_store->set_max_hash_frequency(max_freq);

  // Compute frequency threshold
  std::vector<std::size_t> frequencies;
  frequencies.reserve(result.seed_store->size());
  for (const auto& [hash, hits] : result.seed_store->data()) {
    frequencies.push_back(hits.size());
  }

  std::size_t threshold = max_freq + 1;
  if (!frequencies.empty() && config.seed_freq_cutoff < 1.0) {
    std::sort(frequencies.begin(), frequencies.end());
    double fraction = std::clamp(config.seed_freq_cutoff, 0.0, 1.0);
    std::size_t freq_pos = static_cast<std::size_t>(frequencies.size() * fraction);
    freq_pos = std::min(freq_pos, frequencies.size() - 1);
    threshold = frequencies[freq_pos] + 1;
  }
  result.seed_store->set_frequency_threshold(threshold);

  LOG_INFO("NodeFirstIndex: " + std::to_string(result.seeds_interior) + " interior + " +
           std::to_string(result.seeds_boundary) +
           " boundary = " + std::to_string(result.seeds_interior + result.seeds_boundary) +
           " seeds, " + std::to_string(result.seeds_unique) +
           " unique (max_freq=" + std::to_string(max_freq) + ")");

  return result;
}

}  // namespace piru::index
