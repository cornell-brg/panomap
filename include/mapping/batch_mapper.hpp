/**
 * batch_mapper.hpp
 *
 * BatchMapper config, stats, components, and the mapper class.
 *
 * Related:
 *  - batch_mapper.cpp
 *  - chainer.hpp
 *  - map_result.hpp
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/reads/read_provider.hpp"
#include "io/results/result_writer.hpp"
#include "mapping/chainer.hpp"
#include "mapping/map_result.hpp"
#include "signal/event_pipelines/event_pipeline_factory.hpp"
#include "signal/tokenizers/tokenizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"

namespace piru::mapping {

class SeedLookup {
public:
  SeedLookup(const index::SeedStore* store, std::size_t freq_threshold,
             std::size_t max_total_hits = 0)
      : store_(store), freq_threshold_(freq_threshold), max_total_hits_(max_total_hits) {}

  void lookup(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& out_hits) const;

private:
  const index::SeedStore* store_{nullptr};  // non-owning
  std::size_t freq_threshold_{0};
  std::size_t max_total_hits_{0};  // 0 = unlimited; cap total hits per read
};

struct BatchMapperConfig {
  std::size_t batch_capacity_reads{4000};
  std::size_t batch_capacity_bytes{512 * 1024 * 1024};  // Reserved for future use.
  int num_threads{-1};                                  // -1 = automatic.

  signal::EventPipelineConfig event_pipeline_config{};  // Unified event detection + normalization
  signal::TokenizerConfig tokenizer_config{};
  signal::SeedExtractorConfig seed_config{};
  float diff_filter{0.35f};  // Skip events within diff of last emitted (0 = disabled)
  std::string chainer_backend{"path-chain"};
  cli::Parsed chainer_parsed{};                   // CLI args forwarded to chainer
  const index::SeedStore* seed_store{nullptr};    // non-owning pointer to loaded SeedStore
  const index::GraphStore* graph_store{nullptr};  // non-owning pointer to loaded GraphStore
  /* Linearization coordinates (needed for DP chaining).
   * Non-owning pointer (from in-memory indexing or future deserialization). */
  const std::vector<std::vector<index::LinearCoordinate>>* linearization_coords{nullptr};
  /* Path lengths for anchor bounds checking (parallel to graph paths) */
  const std::vector<std::size_t>* path_lengths{nullptr};
  /* 1D sort coordinates for SortChainer (non-owning, from index) */
  const std::vector<float>* node_1d_coords{nullptr};
  /* Connected component IDs per node (non-owning, from index) */
  const std::vector<std::uint32_t>* component_ids{nullptr};

  /* Result writer for output (non-owning, optional) */
  io::ResultWriter* result_writer{nullptr};

  /* Anchor merging (passed to chainer) */
  bool enable_anchor_merge{true};
  std::size_t pore_k{0};  // Pore model k, passed to chainer for scoring span

  /* Seed lookup limits */
  std::size_t max_total_hits{100000};  // Per-read hit cap (0 = unlimited, default 100k)

  /* Mapping decision.
   * Single-chain reads: anchor count + event/ref ratio gate.
   * Multi-chain reads: weighted standout (RH2-style) vs threshold. */
  int map_min_mapq{2};                // (unused, kept for CLI compatibility)
  float map_w_bestq{0.35f};           // Weight: absolute mapq strength (bestQ/30)
  float map_w_bestmq{0.05f};          // Weight: mapq standout vs mean
  float map_w_bestmc{0.6f};           // Weight: score standout vs mean
  float map_w_threshold{0.45f};       // Weighted sum threshold to accept
  bool no_early_exit{false};          // If true, process all chunks before deciding

  /* Single-chain gate: replaces old fast path.
   * A single chain is accepted if anchors >= min AND event/ref ratio in [lo, hi]. */
  std::size_t map_sc_min_anchors{5};  // Min anchors for single-chain accept
  float map_sc_ratio_lo{0.7f};        // Min event_span/ref_span ratio
  float map_sc_ratio_hi{1.4f};        // Max event_span/ref_span ratio

  /* Fallback: accept strong chains that standout couldn't decide on.
   * Only applies after all chunks exhausted. Uses EMA of accepted chains
   * (mean + variance) to adaptively learn what "good enough" looks like.
   * Threshold = ema_mean - k * ema_std. A chain passes if within k std devs of the mean. */
  double map_fallback_init_score{200.0};      // Initial EMA score (before any data)
  double map_fallback_init_anchors{20.0};     // Initial EMA anchors (before any data)
  float map_fallback_alpha{0.02f};            // EMA smoothing (0.02 ~ 50-read memory)
  float map_fallback_k{1.0f};                 // Threshold = mean - k * std (higher = stricter)
  bool map_fallback_adaptive{true};           // false = use fixed init values only
};

struct BatchMapperStats {
  std::size_t batches{0};
  std::size_t reads_processed{0};
  std::size_t reads_mapped{0};     // reads with at least one chain
  std::size_t reads_unmapped{0};   // reads with no chains
  std::size_t results_written{0};  // total alignment results written
  std::size_t primary_alignments{0};
  std::size_t secondary_alignments{0};
};

struct BatchBuffer {
  std::vector<io::RawRead> raw_reads;
  std::vector<signal::NormalizedSignal> normalized;
  std::vector<signal::TokenizedSignal> tokenized;
  std::vector<signal::SeedBuffer> seeds;
  std::vector<std::vector<NodeAnchor>> seed_hits;
  std::vector<ReadMapResult> map_results;
  std::size_t num_reads{0};

  void resize(std::size_t capacity);
  void clear();
};

struct PipelineComponents {
  signal::EventPipelinePtr event_pipeline;
  signal::TokenizerPtr tokenizer;
  signal::SeedExtractorPtr seed_extractor;
  const index::SeedStore* seed_store{nullptr};    // non-owning; loaded index
  const index::GraphStore* graph_store{nullptr};  // non-owning; loaded index
  SeedLookup lookup{nullptr, 0};
  ChainerPtr chainer;
};

class BatchMapper {
public:
  BatchMapper(io::ReadProvider& provider, BatchMapperConfig config, std::ostream& output);

  BatchMapperStats process_all();

private:
  void load_batch(BatchBuffer& batch);
  void process_batch(BatchBuffer& batch);
  BatchMapperStats output_batch(const BatchBuffer& batch) const;
  void process_read(BatchBuffer& batch, std::size_t index);
  void lookup_seed_hits(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& hits_out) const;
  PipelineComponents create_components() const;

  // Update EMA with a standout-accepted chain at a given chunk depth
  void recordAcceptedChain(double score, std::size_t anchors, std::size_t chunks_processed);
  // Get adaptive fallback thresholds for a given chunk depth
  void getAdaptiveThresholds(std::size_t chunks_processed, double& out_score,
                             std::size_t& out_anchors) const;

  BatchMapperConfig config_;
  io::ReadProvider& provider_;
  std::unique_ptr<concurrency::Executor> executor_;
  PipelineComponents components_;
  std::ostream& output_;

  // Adaptive fallback: per-chunk-depth EMA of accepted chain stats.
  // Index 0 unused (ck=0 means no chunks processed). Indices 1..max_chunks active.
  // Tracks both mean and variance (exponentially-weighted) of score and anchors.
  mutable std::mutex adaptive_mutex_;
  mutable std::vector<double> ema_score_per_ck_;
  mutable std::vector<double> ema_score_var_per_ck_;
  mutable std::vector<double> ema_anchors_per_ck_;
  mutable std::vector<double> ema_anchors_var_per_ck_;
};

}  // namespace piru::mapping
