/**
 * bucket_indexer.cpp
 *
 * Path-walk bucket indexer: walks embedded paths to generate seeds,
 * scatters them into hash-partitioned buckets, finalizes each bucket
 * independently into a BucketSeedStore.
 *
 * Path-walk is the correct semantic model for both k-mer and minimizer
 * extraction: seed selection depends on path context (the full token
 * stream along a path), not just local node content.
 *
 * Build flow:
 *   1. For each path: squigglize full path, quantize, extract seeds,
 *      map positions to (node_id, offset), push to bucket[hash & mask].
 *      Also record linearization coordinates.
 *   2. Per-bucket (parallel): sort by (hash, node_id, offset), dedup
 *      exact duplicates (shared nodes generate identical hits from
 *      multiple paths), finalize into BucketSeedStore.
 *
 * Threading: TBB parallel_for over paths with thread-local bucket
 * buffers, merged after all paths. Per-bucket finalization is parallel.
 *
 * Related:
 *  - bucket_indexer.hpp
 *  - bucket_seed_store.hpp (final index structure)
 *  - path_walk_indexer.cpp (predecessor, same seed generation pattern)
 *
 * SPDX-License-Identifier: MIT
 */

#include "signal/index/bucket_indexer.hpp"

#include "signal/diff_filter.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#ifdef PIRU_USE_TBB
#include <tbb/enumerable_thread_specific.h>
#endif

#include "core/index/bucket_seed_store.hpp"
#include "core/util/logging.hpp"
#include "core/util/timing.hpp"

namespace piru::index {

namespace {

struct SeedEntry16 {
  std::uint64_t hash;
  std::uint32_t node_id;
  std::uint32_t offset;

  bool operator<(const SeedEntry16& o) const {
    if (hash != o.hash) return hash < o.hash;
    if (node_id != o.node_id) return node_id < o.node_id;
    return offset < o.offset;
  }
  bool operator==(const SeedEntry16& o) const {
    return hash == o.hash && node_id == o.node_id && offset == o.offset;
  }
};

using BucketVec = std::vector<SeedEntry16>;

/* Track where each node starts in the concatenated path sequence. */
struct NodeBoundary {
  std::uint32_t node_id;
  std::size_t base_start;
};

/* Map signal position back to (node_id, local_offset). */
std::pair<std::uint32_t, std::uint32_t> signalPosToNodeOffset(
    std::size_t signal_pos, const std::vector<NodeBoundary>& boundaries) {
  auto it =
      std::upper_bound(boundaries.begin(), boundaries.end(), signal_pos,
                       [](std::size_t pos, const NodeBoundary& b) { return pos < b.base_start; });
  if (it != boundaries.begin()) --it;
  return {it->node_id, static_cast<std::uint32_t>(signal_pos - it->base_start)};
}

}  // namespace

BucketIndexResult bucketIndex(const FlatGraph& graph, const io::KmerModel& model,
                              const signal::Tokenizer& tokenizer,
                              const signal::SeedExtractor& extractor,
                              const BucketIndexConfig& config) {
  BucketIndexResult result;
  result.linearization_coords.resize(graph.nodeCount());
  result.path_lengths.resize(graph.pathCount(), 0);

  const int pore_k = model.k();
  const auto pore_flat = model.buildFlatLookup();
  const std::uint64_t kmer_mask = (1ULL << (2 * pore_k)) - 1;

  const std::uint32_t bucket_bits = 14;
  const std::size_t num_buckets = 1ULL << bucket_bits;
  const std::uint64_t bucket_mask = num_buckets - 1;

  std::vector<BucketVec> buckets(num_buckets);
  std::atomic<std::size_t> atomic_seeds_total{0};

  using LinearCoordVec = std::vector<std::vector<LinearCoordinate>>;

  // =========================================================================
  // Phase 1: Walk paths, extract seeds, scatter to buckets
  // =========================================================================
  timing::start("bucket:path_walk");

  auto process_path = [&](std::size_t path_idx, LinearCoordVec& coords_out,
                          std::vector<BucketVec>& local_buckets) {
    const auto* steps = graph.pathStepsBegin(static_cast<std::uint32_t>(path_idx));
    const std::size_t num_steps = graph.pathStepCount(static_cast<std::uint32_t>(path_idx));

    /* Build boundaries and compute path length */
    std::vector<NodeBoundary> boundaries;
    std::size_t path_len = 0;
    for (std::size_t s = 0; s < num_steps; ++s) {
      std::uint32_t nid = steps[s];
      if (nid >= graph.nodeCount()) continue;
      boundaries.push_back({nid, path_len});
      coords_out[nid].emplace_back(path_idx, static_cast<std::int64_t>(path_len));
      path_len += graph.seqLen(nid);
    }

    if (path_len < static_cast<std::size_t>(pore_k)) return;
    // Safe: each path writes its own slot, no concurrent writes to same index.
    result.path_lengths[path_idx] = path_len;

    /* Squigglize full path using rolling 2-bit kmer */
    const std::size_t num_positions = path_len - static_cast<std::size_t>(pore_k) + 1;
    std::vector<float> raw_signal(num_positions);

    std::uint64_t kmer = 0;
    int valid_bases = 0;
    std::size_t sig_idx = 0;
    std::size_t base_pos = 0;

    for (std::size_t s = 0; s < num_steps; ++s) {
      std::uint32_t nid = steps[s];
      if (nid >= graph.nodeCount()) continue;
      std::uint32_t node_len = static_cast<std::uint32_t>(graph.seqLen(nid));

      for (std::uint32_t j = 0; j < node_len; ++j) {
        if (graph.isN(nid, j)) {
          kmer = 0;
          valid_bases = 0;
        } else {
          std::uint8_t base = graph.base2bit(nid, j);
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

    /* Check for valid signal */
    bool has_valid = false;
    for (std::size_t i = 0; i < num_positions; ++i) {
      if (!std::isnan(raw_signal[i])) {
        has_valid = true;
        break;
      }
    }
    if (!has_valid) return;

    /* Diff filter + quantize + extract seeds */
    signal::NormalizedSignal normalized;
    normalized.samples = std::move(raw_signal);
    signal::apply_diff_filter(normalized, config.diff_filter);
    auto tokenized = tokenizer.quantize(normalized);
    auto seeds = extractor.extract(tokenized);

    if (seeds.seeds.empty()) return;

    /* Scatter to buckets */
    atomic_seeds_total.fetch_add(seeds.seeds.size(), std::memory_order_relaxed);

    for (const auto& seed : seeds.seeds) {
      auto [nid, local_off] = signalPosToNodeOffset(seed.position, boundaries);
      std::size_t bi = seed.hash & bucket_mask;
      local_buckets[bi].push_back({seed.hash, nid, local_off});
    }
  };

#ifdef PIRU_USE_TBB
  if (config.executor) {
    using ThreadBuckets = std::vector<BucketVec>;
    tbb::enumerable_thread_specific<LinearCoordVec> thread_coords(
        [&]() { return LinearCoordVec(graph.nodeCount()); });
    tbb::enumerable_thread_specific<ThreadBuckets> thread_buckets(
        [&]() { return ThreadBuckets(num_buckets); });

    config.executor->parallel_for(std::size_t{0}, static_cast<std::size_t>(graph.pathCount()),
                                  std::size_t{1}, [&](std::size_t path_idx) {
                                    process_path(path_idx, thread_coords.local(),
                                                 thread_buckets.local());
                                  });

    for (auto& coords : thread_coords) {
      for (std::size_t i = 0; i < coords.size(); ++i) {
        auto& dst = result.linearization_coords[i];
        dst.insert(dst.end(), coords[i].begin(), coords[i].end());
      }
    }
    for (auto& tb : thread_buckets) {
      for (std::size_t i = 0; i < num_buckets; ++i) {
        if (!tb[i].empty()) {
          buckets[i].insert(buckets[i].end(), tb[i].begin(), tb[i].end());
          BucketVec().swap(tb[i]);
        }
      }
    }
  } else
#endif
  {
    for (std::size_t path_idx = 0; path_idx < graph.pathCount(); ++path_idx) {
      process_path(path_idx, result.linearization_coords, buckets);
    }
  }

  // seeds_interior/boundary are legacy from the node-first design.
  // Path-walk produces all seeds in one pass. Store total in seeds_interior
  // for backward compatibility with callers that read it.
  result.seeds_interior = atomic_seeds_total.load();
  result.seeds_boundary = 0;
  timing::stop("bucket:path_walk");

  // =========================================================================
  // Phase 2: Per-bucket finalization -> BucketSeedStore
  // =========================================================================
  timing::start("bucket:finalize");

  /* Finalize each bucket: sort, dedup, build.
   * No index-time frequency filtering -- that happens at runtime
   * via --seed-freq-cap. Index stores everything after dedup. */
  std::size_t max_freq = 0;
  std::vector<Bucket> finalized(num_buckets);
  std::vector<std::size_t> per_bucket_max_freq(num_buckets, 0);

  auto finalize_one = [&](std::size_t bi) {
    auto& bv = buckets[bi];
    if (bv.empty()) return;

    std::sort(bv.begin(), bv.end());

    // Dedup exact (hash, node_id, offset) duplicates from shared nodes
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

    // Build (hash, SeedEntry) pairs
    std::vector<std::pair<std::uint64_t, SeedEntry>> pairs;
    pairs.reserve(bv.size());
    std::size_t local_max_freq = 0;

    for (std::size_t i = 0; i < bv.size();) {
      std::size_t j = i + 1;
      while (j < bv.size() && bv[j].hash == bv[i].hash) ++j;
      local_max_freq = std::max(local_max_freq, j - i);
      for (std::size_t k = i; k < j; ++k) {
        pairs.push_back({bv[k].hash, SeedEntry{bv[k].node_id, bv[k].offset}});
      }
      i = j;
    }

    BucketVec().swap(bv);
    finalized[bi] = BucketSeedStore::finalize_bucket(pairs);
    per_bucket_max_freq[bi] = local_max_freq;
  };

#ifdef PIRU_USE_TBB
  if (config.executor) {
    config.executor->parallel_for(std::size_t{0}, num_buckets, std::size_t{1}, finalize_one);
  } else
#endif
  {
    for (std::size_t bi = 0; bi < num_buckets; ++bi) {
      finalize_one(bi);
    }
  }

  // Reduce per-bucket max freq
  for (std::size_t bi = 0; bi < num_buckets; ++bi) {
    max_freq = std::max(max_freq, per_bucket_max_freq[bi]);
  }

  // Count totals
  std::size_t total_unique = 0;
  std::size_t total_entries = 0;
  for (const auto& b : finalized) {
    total_unique += b.keys.size();
    total_entries += b.entries.size();
  }

  // Build extractor metadata
  std::map<std::string, std::string> params;
  const auto& cfg = extractor.config();
  params["backend"] = cfg.backend;
  params["k"] = std::to_string(cfg.k);
  params["stride"] = std::to_string(cfg.stride);
  params["qbits"] = std::to_string(cfg.qbits);
  params["window"] = std::to_string(cfg.window);

  result.seed_store =
      std::make_unique<BucketSeedStore>(std::move(finalized), bucket_bits, extractor.name(),
                                        std::move(params), max_freq, max_freq + 1, 1.0);

  result.seeds_total = total_entries;

  timing::stop("bucket:finalize");

  LOG_INFO("BucketIndex: " + std::to_string(atomic_seeds_total.load()) + " raw seeds, " +
           std::to_string(total_unique) + " unique hashes, " + std::to_string(total_entries) +
           " entries after dedup" + " (max_freq=" + std::to_string(max_freq) +
           ", buckets=" + std::to_string(num_buckets) + ")");

  return result;
}

}  // namespace piru::index
