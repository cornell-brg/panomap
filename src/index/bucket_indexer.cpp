/**
 * bucket_indexer.cpp
 *
 * Bucket-partitioned indexer: processes node interiors once (node-first
 * style), fills boundary gaps via path walks, scatters all seeds into
 * hash-partitioned concurrent queues drained by consumer threads, then
 * sorts and builds FlatSeedStore directly. No HashSeedStore, no
 * unordered_map.
 *
 * Threading model:
 *  - N consumer std::threads, each drains its own ConcurrentQueue
 *  - TBB parallel_for for producers (node interiors + path walks)
 *  - Backpressure: atomic counter, producers yield when queues too deep
 *  - Sequential fallback when no TBB / no executor
 *
 * Related:
 *  - bucket_indexer.hpp
 *  - path_walk_indexer.cpp (rolling 2-bit kmer pattern)
 *  - node_first_indexer.cpp (boundary fill logic)
 *
 * SPDX-License-Identifier: MIT
 */

#include "index/bucket_indexer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include "concurrentqueue.h"

#ifdef PIRU_USE_TBB
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_sort.h>
#endif

#include "util/logging.hpp"
#include "util/timing.hpp"

namespace piru::index {

namespace {

/* 16-byte entry scattered to buckets. Sorted by hash in phase 3. */
struct SeedQueueEntry {
    std::uint64_t hash;
    std::uint32_t node_id;
    std::uint32_t offset;

    bool operator<(const SeedQueueEntry& o) const {
        if (hash != o.hash) return hash < o.hash;
        if (node_id != o.node_id) return node_id < o.node_id;
        return offset < o.offset;
    }
};

/**
 * Squigglize a single node using rolling 2-bit kmer + flat pore model.
 * Returns raw signal values (pore model means, no normalization).
 * Empty result if node is too short or has no valid kmers.
 */
std::vector<float> squigglizeNode(const FlatGraph& graph, std::uint32_t node_id,
                                  const std::vector<float>& pore_flat, int pore_k,
                                  std::uint64_t kmer_mask) {
    const std::size_t node_len = graph.seqLen(node_id);
    if (node_len < static_cast<std::size_t>(pore_k)) return {};

    const std::size_t num_positions = node_len - static_cast<std::size_t>(pore_k) + 1;
    std::vector<float> raw_signal(num_positions);

    std::uint64_t kmer = 0;
    int valid_bases = 0;
    std::size_t sig_idx = 0;

    for (std::uint32_t j = 0; j < static_cast<std::uint32_t>(node_len); ++j) {
        if (graph.isN(node_id, j)) {
            kmer = 0;
            valid_bases = 0;
        } else {
            std::uint8_t base = graph.base2bit(node_id, j);
            kmer = ((kmer << 2) | base) & kmer_mask;
            valid_bases++;
        }

        if (j >= static_cast<std::uint32_t>(pore_k) - 1) {
            if (valid_bases >= pore_k) {
                raw_signal[sig_idx] = pore_flat[kmer];
            } else {
                raw_signal[sig_idx] = std::numeric_limits<float>::quiet_NaN();
            }
            sig_idx++;
        }
    }

    bool has_valid = false;
    for (std::size_t i = 0; i < num_positions; ++i) {
        if (!std::isnan(raw_signal[i])) { has_valid = true; break; }
    }
    if (!has_valid) return {};

    return raw_signal;
}

/**
 * Build a boundary context window across node junctions using fast
 * 2-bit rolling kmer. Replaces the string-based getHashWindow().
 */
std::vector<std::int16_t> getHashWindowFast(
    const FlatGraph& graph, const std::vector<float>& pore_flat,
    std::uint64_t kmer_mask, int pore_k,
    const signal::FuzzyQuantizer& fuzzy_quantizer,
    std::size_t path_idx, std::size_t node_step_idx, std::size_t offset,
    std::size_t window_size) {

    const auto* steps = graph.pathStepsBegin(static_cast<std::uint32_t>(path_idx));
    const std::size_t num_steps = graph.pathStepCount(static_cast<std::uint32_t>(path_idx));

    const std::size_t bases_needed = window_size + static_cast<std::size_t>(pore_k) - 1;

    std::vector<std::uint8_t> bases;
    std::vector<bool> is_n;
    bases.reserve(bases_needed);
    is_n.reserve(bases_needed);

    std::size_t cur_step = node_step_idx;
    std::size_t cur_offset = offset;

    while (bases.size() < bases_needed && cur_step < num_steps) {
        std::uint32_t nid = steps[cur_step];
        if (nid >= graph.nodeCount()) {
            ++cur_step;
            cur_offset = 0;
            continue;
        }

        const std::size_t nlen = graph.seqLen(nid);
        std::size_t avail = (cur_offset < nlen) ? (nlen - cur_offset) : 0;
        std::size_t to_take = std::min(avail, bases_needed - bases.size());

        for (std::size_t b = 0; b < to_take; ++b) {
            std::uint32_t pos = static_cast<std::uint32_t>(cur_offset + b);
            is_n.push_back(graph.isN(nid, pos));
            bases.push_back(graph.base2bit(nid, pos));
        }

        ++cur_step;
        cur_offset = 0;
    }

    if (bases.size() < bases_needed) return {};

    std::vector<float> raw_values;
    raw_values.reserve(window_size);

    std::uint64_t kmer = 0;
    int valid_bases = 0;

    for (std::size_t i = 0; i < bases_needed; ++i) {
        if (is_n[i]) {
            kmer = 0;
            valid_bases = 0;
        } else {
            kmer = ((kmer << 2) | bases[i]) & kmer_mask;
            valid_bases++;
        }

        if (i >= static_cast<std::size_t>(pore_k) - 1) {
            if (valid_bases >= pore_k) {
                raw_values.push_back(pore_flat[kmer]);
            } else {
                return {};
            }
        }
    }

    if (raw_values.size() < window_size) return {};

    signal::NormalizedSignal normalized;
    normalized.samples = std::move(raw_values);
    auto fuzzy = fuzzy_quantizer.quantize(normalized);
    return std::move(fuzzy.tokens);
}

/* Consumer thread: drains one ConcurrentQueue into a flat vector. */
void consumeBucket(moodycamel::ConcurrentQueue<SeedQueueEntry>& queue,
                   std::vector<SeedQueueEntry>& bucket,
                   std::atomic<std::size_t>& queued_count,
                   std::atomic<bool>& producers_done) {
    static constexpr std::size_t kBatchSize = 1024;
    SeedQueueEntry batch[kBatchSize];

    while (true) {
        std::size_t count = queue.try_dequeue_bulk(batch, kBatchSize);
        if (count > 0) {
            bucket.insert(bucket.end(), batch, batch + count);
            queued_count.fetch_sub(count, std::memory_order_relaxed);
        } else if (producers_done.load(std::memory_order_acquire)) {
            // Final drain
            while ((count = queue.try_dequeue_bulk(batch, kBatchSize)) > 0) {
                bucket.insert(bucket.end(), batch, batch + count);
                queued_count.fetch_sub(count, std::memory_order_relaxed);
            }
            break;
        } else {
            std::this_thread::yield();
        }
    }
}

}  // namespace

BucketIndexResult bucketIndex(const FlatGraph& graph, const io::KmerModel& model,
                              const signal::FuzzyQuantizer& fuzzy_quantizer,
                              const signal::SeedExtractor& extractor,
                              const BucketIndexConfig& config) {
    BucketIndexResult result;
    result.linearization_coords.resize(graph.nodeCount());
    result.path_lengths.resize(graph.pathCount(), 0);

    const int pore_k = model.k();
    const auto pore_flat = model.buildFlatLookup();
    const std::uint64_t kmer_mask = (1ULL << (2 * pore_k)) - 1;
    const std::size_t seed_k = config.seed_k;

    const std::size_t num_buckets = config.num_buckets > 0
        ? config.num_buckets
        : (config.executor
            ? static_cast<std::size_t>(config.executor->max_concurrency())
            : 1);

    // Bucket storage
    std::vector<std::vector<SeedQueueEntry>> buckets(num_buckets);
    std::size_t est_seeds_per_bucket = graph.totalBases()
        / std::max(config.seed_stride, std::size_t{1}) / num_buckets;
    for (auto& b : buckets) {
        b.reserve(std::min(est_seeds_per_bucket, std::size_t{64 * 1024 * 1024}));
    }

    // Concurrent queues + consumer threads (only when parallel)
    std::vector<moodycamel::ConcurrentQueue<SeedQueueEntry>> queues(num_buckets);
    std::vector<std::thread> consumers;
    std::atomic<std::size_t> queued_count{0};
    std::atomic<bool> producers_done{false};
    const std::size_t max_queued = config.max_queued_entries;

    const bool use_parallel = (config.executor != nullptr);

    if (use_parallel) {
        for (std::size_t i = 0; i < num_buckets; ++i) {
            consumers.emplace_back(consumeBucket,
                std::ref(queues[i]), std::ref(buckets[i]),
                std::ref(queued_count), std::ref(producers_done));
        }
    }

    std::atomic<std::size_t> atomic_seeds_interior{0};

    // =========================================================================
    // Phase 1: Node interiors
    // =========================================================================
    timing::start("bucket:node_interiors");

    /* Scatter seeds from a node to queues (parallel) or directly to buckets (sequential). */
    auto scatter_seeds = [&](const signal::SeedBuffer& seeds, std::uint32_t node_id) {
        for (const auto& seed : seeds.seeds) {
            std::size_t bucket_idx = seed.hash % num_buckets;
            SeedQueueEntry entry{seed.hash, node_id,
                                 static_cast<std::uint32_t>(seed.position)};
            if (use_parallel) {
                queues[bucket_idx].enqueue(entry);
                queued_count.fetch_add(1, std::memory_order_relaxed);
            } else {
                buckets[bucket_idx].push_back(entry);
            }
        }
    };

    auto process_node = [&](std::size_t node_id) {
        auto raw_signal = squigglizeNode(graph, static_cast<std::uint32_t>(node_id),
                                         pore_flat, pore_k, kmer_mask);
        if (raw_signal.empty()) return;

        signal::NormalizedSignal normalized;
        normalized.samples = std::move(raw_signal);
        auto fuzzy = fuzzy_quantizer.quantize(normalized);

        if (fuzzy.tokens.size() < seed_k) return;

        auto seeds = extractor.extract(fuzzy);
        if (seeds.seeds.empty()) return;

        atomic_seeds_interior.fetch_add(seeds.seeds.size(), std::memory_order_relaxed);
        scatter_seeds(seeds, static_cast<std::uint32_t>(node_id));

        // Backpressure: yield if queues are too deep
        if (use_parallel) {
            while (queued_count.load(std::memory_order_relaxed) > max_queued) {
                std::this_thread::yield();
            }
        }
    };

#ifdef PIRU_USE_TBB
    if (config.executor) {
        config.executor->parallel_for(std::size_t{0}, static_cast<std::size_t>(graph.nodeCount()),
                                      std::size_t{100}, process_node);
    } else
#endif
    {
        for (std::size_t node_id = 0; node_id < graph.nodeCount(); ++node_id) {
            process_node(node_id);
        }
    }

    result.seeds_interior = atomic_seeds_interior.load();
    timing::stop("bucket:node_interiors");

    // =========================================================================
    // Phase 2: Boundary fill (path walks) + linearization coords
    // =========================================================================
    timing::start("bucket:boundary_fill");

    std::atomic<std::size_t> atomic_seeds_boundary{0};

    using LinearCoordVec = std::vector<std::vector<LinearCoordinate>>;

    auto process_path = [&](std::size_t path_idx, LinearCoordVec& coords_out) {
        const auto* steps = graph.pathStepsBegin(static_cast<std::uint32_t>(path_idx));
        const std::size_t num_steps = graph.pathStepCount(static_cast<std::uint32_t>(path_idx));

        std::size_t path_base_pos = 0;
        std::size_t local_seed_count = 0;

        for (std::size_t step_idx = 0; step_idx < num_steps; ++step_idx) {
            std::uint32_t node_id = steps[step_idx];
            if (node_id >= graph.nodeCount()) continue;

            const std::size_t node_len = graph.seqLen(node_id);

            // Record linear coordinate
            coords_out[node_id].emplace_back(
                path_idx, static_cast<std::int64_t>(path_base_pos));

            // Compute where interior indexing covered up to
            std::size_t first_boundary_pos = 0;
            if (node_len >= static_cast<std::size_t>(pore_k) + seed_k - 1) {
                std::size_t num_fuzzy = node_len - static_cast<std::size_t>(pore_k) + 1;
                if (num_fuzzy >= seed_k) {
                    first_boundary_pos = num_fuzzy - seed_k + 1;
                }
            }

            // Fill boundary positions
            for (std::size_t pos = first_boundary_pos; pos < node_len;
                 pos += config.seed_stride) {
                auto window = getHashWindowFast(
                    graph, pore_flat, kmer_mask, pore_k, fuzzy_quantizer,
                    path_idx, step_idx, pos, seed_k);

                if (window.empty()) continue;

                signal::FuzzyQuantizedSignal fuzzy_window;
                fuzzy_window.tokens = std::move(window);

                auto seeds = extractor.extract(fuzzy_window);
                for (const auto& seed : seeds.seeds) {
                    std::size_t bucket_idx = seed.hash % num_buckets;
                    SeedQueueEntry entry{seed.hash, node_id,
                                         static_cast<std::uint32_t>(pos + seed.position)};
                    if (use_parallel) {
                        queues[bucket_idx].enqueue(entry);
                        queued_count.fetch_add(1, std::memory_order_relaxed);
                    } else {
                        buckets[bucket_idx].push_back(entry);
                    }
                    ++local_seed_count;
                }
            }

            path_base_pos += node_len;
        }

        result.path_lengths[path_idx] = path_base_pos;

        if (local_seed_count > 0) {
            atomic_seeds_boundary.fetch_add(local_seed_count, std::memory_order_relaxed);
        }
    };

#ifdef PIRU_USE_TBB
    if (config.executor) {
        tbb::enumerable_thread_specific<LinearCoordVec> thread_coords(
            [&]() { return LinearCoordVec(graph.nodeCount()); });

        config.executor->parallel_for(std::size_t{0}, static_cast<std::size_t>(graph.pathCount()),
                                      std::size_t{1}, [&](std::size_t path_idx) {
            process_path(path_idx, thread_coords.local());
        });

        // Merge thread-local linearization coords
        for (auto& coords : thread_coords) {
            for (std::size_t i = 0; i < coords.size(); ++i) {
                auto& dst = result.linearization_coords[i];
                dst.insert(dst.end(), coords[i].begin(), coords[i].end());
            }
        }
    } else
#endif
    {
        for (std::size_t path_idx = 0; path_idx < graph.pathCount(); ++path_idx) {
            process_path(path_idx, result.linearization_coords);
        }
    }

    result.seeds_boundary = atomic_seeds_boundary.load();
    timing::stop("bucket:boundary_fill");

    // Signal consumers to finish and wait
    if (use_parallel) {
        producers_done.store(true, std::memory_order_release);
        for (auto& t : consumers) t.join();
    }

    // =========================================================================
    // Phase 3: Sort + build FlatSeedStore
    // =========================================================================
    timing::start("bucket:finalize");

    // Concatenate all buckets into one flat array
    std::size_t total_entries = 0;
    for (const auto& b : buckets) total_entries += b.size();

    std::vector<SeedQueueEntry> all_entries;
    all_entries.reserve(total_entries);
    for (auto& b : buckets) {
        all_entries.insert(all_entries.end(), b.begin(), b.end());
        std::vector<SeedQueueEntry>().swap(b);
    }
    buckets.clear();

    // Sort by hash (then node_id, offset for determinism)
#ifdef PIRU_USE_TBB
    if (config.executor) {
        tbb::parallel_sort(all_entries.begin(), all_entries.end());
    } else
#endif
    {
        std::sort(all_entries.begin(), all_entries.end());
    }

    // Global frequency filter
    std::size_t freq_threshold = std::numeric_limits<std::size_t>::max();
    if (config.seed_freq_cutoff < 1.0 && !all_entries.empty()) {
        std::vector<std::size_t> freqs;
        std::size_t run_start = 0;
        for (std::size_t i = 1; i <= all_entries.size(); ++i) {
            if (i == all_entries.size() || all_entries[i].hash != all_entries[run_start].hash) {
                freqs.push_back(i - run_start);
                run_start = i;
            }
        }
        std::sort(freqs.begin(), freqs.end());
        double frac = std::clamp(config.seed_freq_cutoff, 0.0, 1.0);
        std::size_t pos = std::min(
            static_cast<std::size_t>(freqs.size() * frac), freqs.size() - 1);
        freq_threshold = freqs[pos] + 1;
    }

    // Build CSR arrays
    std::vector<std::uint64_t> csr_hashes;
    std::vector<std::uint32_t> csr_offsets;
    std::vector<SeedEntry> csr_entries;

    csr_hashes.reserve(total_entries / 4);
    csr_offsets.reserve(total_entries / 4 + 1);
    csr_entries.reserve(total_entries);

    std::size_t max_freq = 0;
    std::size_t run_start = 0;

    for (std::size_t i = 0; i <= all_entries.size(); ++i) {
        bool end_of_run = (i == all_entries.size()) ||
                          (i > 0 && all_entries[i].hash != all_entries[run_start].hash);

        if (end_of_run && i > 0) {
            std::size_t run_len = i - run_start;
            std::uint64_t hash = all_entries[run_start].hash;

            if (run_len < freq_threshold) {
                csr_hashes.push_back(hash);
                csr_offsets.push_back(static_cast<std::uint32_t>(csr_entries.size()));
                for (std::size_t j = run_start; j < i; ++j) {
                    csr_entries.push_back(SeedEntry{
                        all_entries[j].node_id, all_entries[j].offset});
                }
                max_freq = std::max(max_freq, run_len);
            }

            run_start = i;
        }
    }

    csr_offsets.push_back(static_cast<std::uint32_t>(csr_entries.size()));

    all_entries.clear();
    all_entries.shrink_to_fit();

    result.seeds_total = csr_entries.size();

    // Build extractor metadata
    std::map<std::string, std::string> params;
    const auto& cfg = extractor.config();
    params["backend"] = cfg.backend;
    params["k"] = std::to_string(cfg.k);
    params["stride"] = std::to_string(cfg.stride);
    params["qbits"] = std::to_string(cfg.qbits);
    params["window"] = std::to_string(cfg.window);

    result.seed_store = std::make_unique<FlatSeedStore>(
        std::move(csr_hashes), std::move(csr_offsets), std::move(csr_entries),
        extractor.name(), std::move(params), max_freq, freq_threshold,
        config.seed_freq_cutoff);

    timing::stop("bucket:finalize");

    LOG_INFO("BucketIndex: " + std::to_string(result.seeds_interior) + " interior + " +
             std::to_string(result.seeds_boundary) + " boundary = " +
             std::to_string(result.seeds_interior + result.seeds_boundary) +
             " total, " + std::to_string(result.seed_store->size()) +
             " unique hashes (max_freq=" + std::to_string(max_freq) +
             ", freq_threshold=" + std::to_string(freq_threshold) +
             ", buckets=" + std::to_string(num_buckets) + ")");

    return result;
}

}  // namespace piru::index
