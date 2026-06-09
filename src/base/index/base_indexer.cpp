/**
 * base_indexer.cpp
 *
 * Implementation of the base-mode path-walk minimizer indexer.
 *
 * Algorithm (per path, parallel):
 *   1. Walk path steps; record (node_id, base_start) boundaries and
 *      per-node linearization coords.
 *   2. Roll a 2-bit k-mer over the concatenated path bases; reset on
 *      any N base. Hash via minimap2's hash64. Skip any k-mer whose
 *      window has fewer than k valid bases.
 *   3. Slide a window of size w over the per-position hashes; emit the
 *      leftmost-minimum hash as a seed (deduped against the previous
 *      emitted index). Map the seed start position back to (node_id,
 *      local_offset) via the boundary table; scatter to bucket.
 *
 * Phase 2 (per-bucket, parallel): sort, dedup exact (hash, node, offset)
 * triples that arise from shared nodes, finalize into a BucketSeedStore.
 *
 * Hashing is directional (no canonicalization): the directional FlatGraph
 * already carries forward + reverse-complement copies of each segment, so
 * each strand emits its own minimizer hashes. Reads will hash the same way.
 *
 * Related:
 *  - base_indexer.hpp
 *  - signal/index/bucket_indexer.cpp (parallel signal-mode flow)
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/index/base_indexer.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#ifdef PANOMAP_USE_TBB
#include <tbb/enumerable_thread_specific.h>
#endif

#include "core/index/bucket_seed_store.hpp"
#include "core/util/kmer_hash.hpp"
#include "core/util/logging.hpp"
#include "core/util/timing.hpp"

namespace panomap::base {

namespace {

constexpr std::uint32_t kBucketBits = 14;
constexpr std::uint64_t kInvalidHash = std::numeric_limits<std::uint64_t>::max();

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

struct NodeBoundary {
  std::uint32_t node_id;
  std::size_t base_start;
};

std::pair<std::uint32_t, std::uint32_t> basePosToNodeOffset(
    std::size_t base_pos, const std::vector<NodeBoundary>& boundaries) {
  auto it =
      std::upper_bound(boundaries.begin(), boundaries.end(), base_pos,
                       [](std::size_t pos, const NodeBoundary& b) { return pos < b.base_start; });
  if (it != boundaries.begin()) --it;
  return {it->node_id, static_cast<std::uint32_t>(base_pos - it->base_start)};
}

}  // namespace

BaseBucketIndexResult bucketIndexBase(const panomap::index::FlatGraph& graph,
                                      const BaseBucketIndexConfig& config) {
  using panomap::index::Bucket;
  using panomap::index::BucketSeedStore;
  using panomap::index::LinearCoordinate;
  using panomap::index::SeedEntry;

  BaseBucketIndexResult result;
  result.linearization_coords.resize(graph.nodeCount());
  result.path_lengths.resize(graph.pathCount(), 0);

  const std::size_t k = config.k;
  const std::size_t w = std::max<std::size_t>(config.w, 1);
  if (k == 0 || k > 32) {
    throw std::runtime_error("base_indexer: k must be in [1, 32]");
  }
  const std::uint64_t kmer_mask = (k == 32) ? ~std::uint64_t{0} : ((1ULL << (2 * k)) - 1);

  const std::size_t num_buckets = 1ULL << kBucketBits;
  const std::uint64_t bucket_mask = num_buckets - 1;

  std::vector<BucketVec> buckets(num_buckets);
  std::atomic<std::size_t> atomic_seeds_total{0};

  using LinearCoordVec = std::vector<std::vector<LinearCoordinate>>;

  /* Phase 1: walk paths, extract minimizers, scatter to buckets. */

  timing::start("base_bucket:path_walk");

  auto process_path = [&](std::size_t path_idx, LinearCoordVec& coords_out,
                          std::vector<BucketVec>& local_buckets) {
    const auto* steps = graph.pathStepsBegin(static_cast<std::uint32_t>(path_idx));
    const std::size_t num_steps = graph.pathStepCount(static_cast<std::uint32_t>(path_idx));

    /* Build boundaries + linearization coords; compute path length. */
    std::vector<NodeBoundary> boundaries;
    boundaries.reserve(num_steps);
    std::size_t path_len = 0;
    for (std::size_t s = 0; s < num_steps; ++s) {
      std::uint32_t nid = steps[s];
      if (nid >= graph.nodeCount()) continue;
      boundaries.push_back({nid, path_len});
      coords_out[nid].emplace_back(path_idx, static_cast<std::int64_t>(path_len));
      path_len += graph.seqLen(nid);
    }

    if (path_len < k) return;
    result.path_lengths[path_idx] = path_len;

    /* Roll a directional 2-bit k-mer; emit hash64 per valid k-mer position. */
    const std::size_t num_kmers = path_len - k + 1;
    std::vector<std::uint64_t> hashes(num_kmers, kInvalidHash);

    std::uint64_t kmer = 0;
    int valid_bases = 0;
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
          if (valid_bases < static_cast<int>(k)) ++valid_bases;
        }

        if (base_pos + 1 >= k) {
          const std::size_t kmer_idx = base_pos + 1 - k;
          if (valid_bases >= static_cast<int>(k)) {
            hashes[kmer_idx] = panomap::core::hash64(kmer, kmer_mask);
          }
        }
        ++base_pos;
      }
    }

    /* Sliding window of size w; emit leftmost-min as the minimizer. */
    if (num_kmers < w) return;

    std::size_t last_emitted = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i + w <= num_kmers; ++i) {
      std::size_t min_idx = i;
      std::uint64_t min_hash = hashes[i];
      for (std::size_t j = i + 1; j < i + w; ++j) {
        if (hashes[j] < min_hash) {
          min_hash = hashes[j];
          min_idx = j;
        }
      }
      if (min_hash == kInvalidHash) continue;
      if (min_idx == last_emitted) continue;
      last_emitted = min_idx;

      auto [nid, local_off] = basePosToNodeOffset(min_idx, boundaries);
      std::size_t bi = min_hash & bucket_mask;
      local_buckets[bi].push_back({min_hash, nid, local_off});
      atomic_seeds_total.fetch_add(1, std::memory_order_relaxed);
    }
  };

#ifdef PANOMAP_USE_TBB
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

  timing::stop("base_bucket:path_walk");

  /* Phase 2: per-bucket finalization (sort, dedup, build). */

  timing::start("base_bucket:finalize");

  std::size_t max_freq = 0;
  std::vector<Bucket> finalized(num_buckets);
  std::vector<std::size_t> per_bucket_max_freq(num_buckets, 0);

  auto finalize_one = [&](std::size_t bi) {
    auto& bv = buckets[bi];
    if (bv.empty()) return;

    std::sort(bv.begin(), bv.end());
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

    std::vector<std::pair<std::uint64_t, SeedEntry>> pairs;
    pairs.reserve(bv.size());
    std::size_t local_max_freq = 0;

    for (std::size_t i = 0; i < bv.size();) {
      std::size_t j = i + 1;
      while (j < bv.size() && bv[j].hash == bv[i].hash) ++j;
      local_max_freq = std::max(local_max_freq, j - i);
      for (std::size_t kk = i; kk < j; ++kk) {
        pairs.push_back({bv[kk].hash, SeedEntry{bv[kk].node_id, bv[kk].offset}});
      }
      i = j;
    }

    BucketVec().swap(bv);
    finalized[bi] = BucketSeedStore::finalize_bucket(pairs);
    per_bucket_max_freq[bi] = local_max_freq;
  };

#ifdef PANOMAP_USE_TBB
  if (config.executor) {
    config.executor->parallel_for(std::size_t{0}, num_buckets, std::size_t{1}, finalize_one);
  } else
#endif
  {
    for (std::size_t bi = 0; bi < num_buckets; ++bi) {
      finalize_one(bi);
    }
  }

  for (std::size_t bi = 0; bi < num_buckets; ++bi) {
    max_freq = std::max(max_freq, per_bucket_max_freq[bi]);
  }

  std::size_t total_unique = 0;
  std::size_t total_entries = 0;
  for (const auto& b : finalized) {
    total_unique += b.keys.size();
    total_entries += b.entries.size();
  }

  std::map<std::string, std::string> params;
  params["backend"] = "base_minimizer";
  params["k"] = std::to_string(k);
  params["window"] = std::to_string(w);

  result.seed_store =
      std::make_unique<BucketSeedStore>(std::move(finalized), kBucketBits, "base_minimizer",
                                        std::move(params), max_freq, max_freq + 1, 1.0);

  result.seeds_total = total_entries;

  timing::stop("base_bucket:finalize");

  LOG_INFO("BaseBucketIndex: " + std::to_string(atomic_seeds_total.load()) + " raw minimizers, " +
           std::to_string(total_unique) + " unique hashes, " + std::to_string(total_entries) +
           " entries after dedup (max_freq=" + std::to_string(max_freq) +
           ", buckets=" + std::to_string(num_buckets) + ", k=" + std::to_string(k) +
           ", w=" + std::to_string(w) + ")");

  return result;
}

}  // namespace panomap::base
