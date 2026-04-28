/**
 * base_mapper.hpp
 *
 * Base-mode read mapper. Parallel to signal::mapping::BatchMapper but
 * with the signal-specific stages stripped: per-read flow is
 * (basecalled_sequence -> minimizer seeds -> lookup -> chainer ->
 * decision -> output). No chunking in this first cut -- the whole read
 * is processed at once. Chunked / streaming evaluation is a future
 * concern that mirrors signal mode's chunk loop.
 *
 * Reuses from core: NodeAnchor / Chain / ChainResult, ChainerPtr,
 * ResultWriter, ReadMapResult.
 *
 * Copies from signal/mapping/batch_mapper.cpp (per dev-108 architecture
 * decision: "two mappers share NOTHING at the read-processing layer"):
 *  - SeedLookup -> BaseSeedLookup (takes base::SeedBuffer)
 *  - computeMapq / computeChainSpans / checkMappingDecision (helpers)
 *  - BatchBuffer / process_batch / parallel scaffolding
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "base/io/reads/fastq_provider.hpp"
#include "base/seeder/base_seeder.hpp"
#include "core/concurrency/executor.hpp"
#include "core/index/graph_store.hpp"
#include "core/index/linearizer.hpp"
#include "core/index/seed_store.hpp"
#include "core/io/results/result_writer.hpp"
#include "core/mapping/chainer.hpp"
#include "core/mapping/map_result.hpp"
#include "cli/parse.hpp"

namespace piru::base::mapping {

class BaseSeedLookup {
 public:
  BaseSeedLookup(const piru::index::SeedStore* store, std::size_t mid_occ,
                 std::size_t max_max_occ, std::size_t occ_dist,
                 std::size_t max_total_hits = 0)
      : store_(store),
        mid_occ_(mid_occ),
        max_max_occ_(max_max_occ),
        occ_dist_(occ_dist),
        max_total_hits_(max_total_hits) {}

  // Look up seeds, applying minimap2-style adaptive seed selection:
  //  1. Soft cap: seeds with hits > mid_occ are "high-occ".
  //  2. In each query window of length `occ_dist`, keep up to (window_len /
  //     occ_dist) of the lowest-frequency high-occ seeds; drop the rest.
  //  3. Hard cap: anything with hits > max_max_occ is dropped no matter what.
  // See related/minimap2/seed.c:mm_seed_select.
  void lookup(const piru::base::SeedBuffer& seeds,
              std::vector<piru::mapping::NodeAnchor>& out_hits) const;

 private:
  const piru::index::SeedStore* store_{nullptr};
  std::size_t mid_occ_{0};       // soft cap (high-occ boundary)
  std::size_t max_max_occ_{0};   // hard cap (always-drop)
  std::size_t occ_dist_{0};      // query window for adaptive keep
  std::size_t max_total_hits_{0};
};

struct BaseMapperConfig {
  std::size_t batch_capacity_reads{1000};
  int num_threads{-1};

  /* Seeding */
  piru::base::BaseSeederConfig seeder{};

  /* Chaining */
  std::string chainer_backend{"path-chain"};
  piru::cli::Parsed chainer_parsed{};
  bool enable_anchor_merge{true};

  /* Index views (non-owning) */
  const piru::index::SeedStore* seed_store{nullptr};
  const piru::index::GraphStore* graph_store{nullptr};
  const std::vector<std::vector<piru::index::LinearCoordinate>>* linearization_coords{nullptr};
  const std::vector<std::size_t>* path_lengths{nullptr};
  const std::vector<float>* node_1d_coords{nullptr};
  const std::vector<std::uint32_t>* component_ids{nullptr};

  /* Output */
  piru::io::ResultWriter* result_writer{nullptr};
  std::unordered_map<std::string, piru::io::ResultWriter*> per_file_writers;

  /* Seed-frequency filtering (minimap2-style mid_occ + adaptive in-window
   * selection + max_max_occ hard cap). Defaults match minimap2:
   *  - mid_occ_frac = 2e-4 (top 0.02% of seeds are "high-occ")
   *  - min_mid_occ = 10, max_mid_occ = 1e6 (clamps for the auto formula)
   *  - max_max_occ = 4095 (absolute hard cap)
   *  - occ_dist = 500 (keep ~1 high-occ seed per 500 bp query window)
   *
   * mid_occ < 0 means "auto from frac"; >= 0 overrides. */
  float mid_occ_frac{2e-4f};
  std::size_t min_mid_occ{10};
  std::size_t max_mid_occ{1000000};
  std::int64_t mid_occ_override{-1};
  std::size_t max_max_occ{4095};
  std::size_t occ_dist{500};

  /* Lookup limits */
  std::size_t max_total_hits{100000};

  /* Chunked evaluation (AS-style early decision).
   * chunk_bp = 0 disables chunking; whole read processed at once.
   * Default 450 bp ~= 1 second of R9.4 translocation, matching the
   * signal-mode default of 4000 samples @ 4 kHz. */
  std::size_t chunk_bp{0};
  std::size_t max_chunks{0};   // 0 = unlimited
  bool no_early_exit{false};

  /* Decision policy (mirrors signal mode shape; defaults may need
   * retuning per dev-108 G4). */
  float map_w_bestq{0.35f};
  float map_w_bestmq{0.05f};
  float map_w_bestmc{0.6f};
  float map_w_threshold{0.45f};
  std::size_t map_sc_min_anchors{5};
  float map_sc_ratio_lo{0.7f};
  float map_sc_ratio_hi{1.4f};
  double map_fallback_floor_score{100.0};
  float map_fallback_alpha{0.02f};
  bool map_fallback_adaptive{true};
};

struct BaseMapperStats {
  std::size_t batches{0};
  std::size_t reads_processed{0};
  std::size_t reads_mapped{0};
  std::size_t reads_unmapped{0};
  std::size_t results_written{0};
  std::size_t primary_alignments{0};
  std::size_t secondary_alignments{0};
};

struct BaseBatchBuffer {
  std::vector<piru::base::io::FastqRead> reads;
  std::vector<piru::base::SeedBuffer> seeds;
  std::vector<std::vector<piru::mapping::NodeAnchor>> seed_hits;
  std::vector<piru::mapping::ReadMapResult> map_results;
  std::size_t num_reads{0};

  void resize(std::size_t capacity);
  void clear();
};

struct BasePipelineComponents {
  const piru::index::SeedStore* seed_store{nullptr};
  const piru::index::GraphStore* graph_store{nullptr};
  BaseSeedLookup lookup{nullptr, 0, 0, 0};
  piru::mapping::ChainerPtr chainer;
};

class BaseMapper {
 public:
  BaseMapper(piru::base::io::FastqProvider& provider, BaseMapperConfig config,
             std::ostream& output);

  BaseMapperStats process_all();

 private:
  void load_batch(BaseBatchBuffer& batch);
  void process_batch(BaseBatchBuffer& batch);
  BaseMapperStats output_batch(const BaseBatchBuffer& batch) const;
  void process_read(BaseBatchBuffer& batch, std::size_t index);
  BasePipelineComponents create_components() const;

  // Per-chunk-depth EMA of multi-chain accepted scores. Mirrors signal
  // mode's adaptive fallback: threshold = max(ema_at_ck, floor).
  void recordAcceptedChain(double score, std::size_t chunks_processed);
  double getFallbackThreshold(std::size_t chunks_processed) const;

  BaseMapperConfig config_;
  piru::base::io::FastqProvider& provider_;
  std::unique_ptr<piru::concurrency::Executor> executor_;
  BasePipelineComponents components_;
  std::ostream& output_;

  mutable std::mutex adaptive_mutex_;
  mutable std::vector<double> ema_score_per_ck_;
};

}  // namespace piru::base::mapping
