/**
 * base_mapper.cpp
 *
 * Per-read flow:
 *   1. Read FASTQ record (id + bases) from provider.
 *   2. base_seeder::extract_minimizers -> base::SeedBuffer.
 *   3. BaseSeedLookup -> vector<NodeAnchor>.
 *   4. Chainer::chain -> ChainResult.
 *   5. computeMapq + checkMappingDecision -> ReadMapResult{is_mapped, mappings}.
 *   6. ResultWriter::write.
 *
 * Helpers (computeMapq, computeChainSpans, checkMappingDecision) are
 * copy-pasted from signal/mapping/batch_mapper.cpp -- they are
 * modality-agnostic in shape but were embedded there. dev-108 locks
 * the architecture as "no shared post-seed pipeline"; we honour that
 * by duplicating rather than extracting to core.
 *
 * SPDX-License-Identifier: MIT
 */

#include "base/mapping/base_mapper.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "core/mapping/pan_chainer.hpp"
#include "core/mapping/path_chainer.hpp"
#include "core/mapping/sort_chainer.hpp"
#include "core/util/logging.hpp"

namespace piru::base::mapping {

using piru::mapping::Chain;
using piru::mapping::ChainResult;
using piru::mapping::DecisionPath;
using piru::mapping::Mapping;
using piru::mapping::NodeAnchor;
using piru::mapping::PanChainer;
using piru::mapping::PanChainerConfig;
using piru::mapping::PathChainer;
using piru::mapping::PathChainerConfig;
using piru::mapping::ReadMapResult;
using piru::mapping::SortChainer;
using piru::mapping::SortChainerConfig;

namespace {

// RH2-style mapq. Copy of the signal-side helper.
int computeMapq(double pri_score, double sec_score, std::size_t anchors, int min_chain_sc) {
  if (pri_score <= 0.0) return 0;
  constexpr float q_coef = 40.0f;

  float pen_s1 = static_cast<float>((pri_score > 100 ? 1.0 : 0.01 * pri_score));
  float pen_cm = static_cast<float>(anchors > 10 ? 1.0 : 0.1 * anchors);
  pen_cm = std::min(pen_s1, pen_cm);

  double subsc = std::max(sec_score, static_cast<double>(min_chain_sc));
  float x = static_cast<float>(subsc / pri_score);

  int mapq = static_cast<int>(pen_cm * q_coef * (1.0f - x) * std::log(pri_score));
  mapq = std::max(mapq, 0);
  return std::min(mapq, 60);
}

void computeChainSpans(const Chain& chain, std::uint32_t& event_span, std::int64_t& ref_span) {
  if (chain.anchors.empty()) {
    event_span = 0;
    ref_span = 0;
    return;
  }
  std::uint32_t min_q = chain.anchors.front().read_pos;
  std::uint32_t max_q = chain.anchors.front().read_pos;
  std::int64_t min_r = chain.anchors.front().ref_coord;
  std::int64_t max_r = chain.anchors.front().ref_coord;
  for (const auto& a : chain.anchors) {
    min_q = std::min(min_q, a.read_pos);
    max_q = std::max(max_q, a.read_pos + a.length);
    min_r = std::min(min_r, a.ref_coord);
    max_r = std::max(max_r, a.ref_coord + static_cast<std::int64_t>(a.length));
  }
  event_span = max_q - min_q;
  ref_span = max_r - min_r;
}

bool checkMappingDecision(const std::vector<Chain>& chains, float w_bestq, float w_bestmq,
                          float w_bestmc, float w_threshold, std::size_t sc_min_anchors,
                          float sc_ratio_lo, float sc_ratio_hi, float* out_standout = nullptr,
                          DecisionPath* out_path = nullptr) {
  if (out_path) *out_path = DecisionPath::kUnmapped;
  if (chains.empty()) {
    if (out_standout) *out_standout = 0.0f;
    return false;
  }

  const auto& best = chains[0];
  double sec_score = (chains.size() > 1) ? chains[1].score : 0.0;

  int best_mapq = computeMapq(best.score, sec_score, best.anchors.size(), 15);
  float bestQ = static_cast<float>(best_mapq);
  float bestC = static_cast<float>(best.score);
  float meanQ = 0.0f, meanC = 0.0f;
  for (std::size_t i = 0; i < chains.size(); ++i) {
    double sec = (i + 1 < chains.size()) ? chains[i + 1].score : 0.0;
    int mq = computeMapq(chains[i].score, sec, chains[i].anchors.size(), 15);
    meanQ += static_cast<float>(mq);
    meanC += static_cast<float>(chains[i].score);
  }
  meanQ /= static_cast<float>(chains.size());
  meanC /= static_cast<float>(chains.size());

  float r_bestq = std::min(bestQ / 30.0f, 1.0f);
  float r_bestmq = (bestQ > 0) ? std::max(0.0f, 1.0f - meanQ / bestQ) : 0.0f;
  float r_bestmc = (bestC > 0) ? std::max(0.0f, 1.0f - meanC / bestC) : 0.0f;
  float weighted = w_bestq * r_bestq + w_bestmq * r_bestmq + w_bestmc * r_bestmc;
  if (out_standout) *out_standout = weighted;

  if (chains.size() == 1) {
    if (best.anchors.size() < sc_min_anchors) return false;
    std::uint32_t event_span = 0;
    std::int64_t ref_span = 0;
    computeChainSpans(best, event_span, ref_span);
    if (ref_span <= 0) return false;
    float ratio = static_cast<float>(event_span) / static_cast<float>(ref_span);
    if (ratio >= sc_ratio_lo && ratio <= sc_ratio_hi) {
      if (out_path) *out_path = DecisionPath::kGate;
      return true;
    }
    return false;
  }

  if (weighted >= w_threshold) {
    if (out_path) *out_path = DecisionPath::kStandout;
    return true;
  }
  return false;
}

}  // namespace

void BaseSeedLookup::lookup(const piru::base::SeedBuffer& seeds,
                            std::vector<NodeAnchor>& out_hits) const {
  if (!store_) return;
  out_hits.clear();
  out_hits.reserve(seeds.seeds.size());
  for (const auto& seed : seeds.seeds) {
    auto hits = store_->lookup(seed.hash);
    if (hits.count == 0) continue;
    if (hits.count > freq_threshold_) continue;
    for (std::size_t j = 0; j < hits.count; ++j) {
      const auto& h = hits.data[j];
      out_hits.push_back(NodeAnchor{
          .node_id = h.node_id,
          .offset = h.offset,
          .read_pos = static_cast<std::uint32_t>(seed.position),
          .span = static_cast<std::uint16_t>(std::min(seed.length, std::size_t{0xFFFF})),
          .length = static_cast<std::uint16_t>(std::min(seed.length, std::size_t{0xFFFF})),
      });
    }
    if (max_total_hits_ > 0 && out_hits.size() >= max_total_hits_) break;
  }
}

void BaseBatchBuffer::resize(std::size_t capacity) {
  reads.resize(capacity);
  seeds.resize(capacity);
  seed_hits.resize(capacity);
  map_results.resize(capacity);
  num_reads = 0;
}

void BaseBatchBuffer::clear() {
  for (std::size_t i = 0; i < num_reads; ++i) {
    reads[i] = piru::base::io::FastqRead{};
    seeds[i].seeds.clear();
    seed_hits[i].clear();
    map_results[i] = ReadMapResult{};
  }
  num_reads = 0;
}

BaseMapper::BaseMapper(piru::base::io::FastqProvider& provider, BaseMapperConfig config,
                       std::ostream& output)
    : config_(std::move(config)),
      provider_(provider),
      executor_(piru::concurrency::make_executor(config_.num_threads)),
      components_(create_components()),
      output_(output) {
  // Pre-size per-chunk-depth EMA. Index 0 unused (ck=0 means no chunks
  // processed). Initialised to floor so threshold = max(ema, floor)
  // begins at floor.
  const std::size_t mc = config_.max_chunks;
  const std::size_t sz = (mc > 0 ? mc : 20) + 1;
  ema_score_per_ck_.assign(sz, config_.map_fallback_floor_score);
}

void BaseMapper::recordAcceptedChain(double score, std::size_t chunks_processed) {
  if (!config_.map_fallback_adaptive) return;
  if (chunks_processed == 0 || chunks_processed >= ema_score_per_ck_.size()) return;
  std::lock_guard<std::mutex> lock(adaptive_mutex_);
  const double a = config_.map_fallback_alpha;
  ema_score_per_ck_[chunks_processed] =
      a * score + (1.0 - a) * ema_score_per_ck_[chunks_processed];
}

double BaseMapper::getFallbackThreshold(std::size_t chunks_processed) const {
  const double floor = config_.map_fallback_floor_score;
  if (!config_.map_fallback_adaptive) return floor;
  if (chunks_processed == 0 || chunks_processed >= ema_score_per_ck_.size()) return floor;
  std::lock_guard<std::mutex> lock(adaptive_mutex_);
  return std::max(ema_score_per_ck_[chunks_processed], floor);
}

BasePipelineComponents BaseMapper::create_components() const {
  BasePipelineComponents comps;
  comps.seed_store = config_.seed_store;
  comps.graph_store = config_.graph_store;
  if (!comps.seed_store) {
    throw std::runtime_error("BaseMapper requires a SeedStore for lookup");
  }

  if (config_.chainer_backend == "path-chain") {
    if (!config_.linearization_coords) {
      throw std::runtime_error("PathChainer requires linearization_coords");
    }
    if (!config_.path_lengths) {
      throw std::runtime_error("PathChainer requires path_lengths");
    }
    auto pc = PathChainerConfig::from_parsed(config_.chainer_parsed);
    pc.merge_anchors = config_.enable_anchor_merge;
    pc.pore_k = 0;  // base mode: anchor span = seed_k, not derived from pore model
    comps.chainer = std::make_unique<PathChainer>(pc, *config_.linearization_coords,
                                                  *config_.path_lengths, config_.node_1d_coords,
                                                  config_.component_ids);
  } else if (config_.chainer_backend == "sort-chain") {
    if (!config_.node_1d_coords) {
      throw std::runtime_error("SortChainer requires node_1d_coords");
    }
    auto sc = SortChainerConfig::from_parsed(config_.chainer_parsed);
    sc.pore_k = 0;
    std::vector<std::uint32_t> bp_lens(config_.graph_store->nodeCount());
    for (std::size_t i = 0; i < bp_lens.size(); ++i) {
      bp_lens[i] = static_cast<std::uint32_t>(config_.graph_store->sequenceLen(i));
    }
    static const std::vector<std::uint32_t> kEmptyCC;
    const auto& cc_ids = config_.component_ids ? *config_.component_ids : kEmptyCC;
    comps.chainer =
        std::make_unique<SortChainer>(sc, *config_.node_1d_coords, std::move(bp_lens), cc_ids);
  } else if (config_.chainer_backend == "pan-chain") {
    if (!config_.node_1d_coords || !config_.linearization_coords || !config_.path_lengths) {
      throw std::runtime_error("PanChainer requires node_1d_coords + linearization + path_lengths");
    }
    auto pc = PanChainerConfig::from_parsed(config_.chainer_parsed);
    pc.pore_k = 0;
    comps.chainer = std::make_unique<PanChainer>(pc, *config_.node_1d_coords,
                                                 *config_.linearization_coords,
                                                 *config_.path_lengths);
  } else {
    throw std::runtime_error("Unknown chainer backend: " + config_.chainer_backend);
  }

  const std::size_t freq_threshold = comps.seed_store->frequency_threshold();
  comps.lookup = BaseSeedLookup(comps.seed_store, freq_threshold, config_.max_total_hits);

  LOG_DEBUG("Pipeline: " + comps.chainer->name() + " chaining (base mode)");
  return comps;
}

BaseMapperStats BaseMapper::process_all() {
  BaseMapperStats stats;
  BaseBatchBuffer batch;
  batch.resize(config_.batch_capacity_reads);

  LOG_DEBUG("BaseMapper starting: batch_capacity_reads=" +
            std::to_string(config_.batch_capacity_reads) +
            ", threads=" + std::to_string(executor_->max_concurrency()));

  while (true) {
    load_batch(batch);
    if (batch.num_reads == 0) break;

    ++stats.batches;
    stats.reads_processed += batch.num_reads;
    process_batch(batch);

    auto bs = output_batch(batch);
    stats.reads_mapped += bs.reads_mapped;
    stats.reads_unmapped += bs.reads_unmapped;
    stats.results_written += bs.results_written;
    stats.primary_alignments += bs.primary_alignments;
    stats.secondary_alignments += bs.secondary_alignments;

    batch.clear();
  }

  LOG_INFO("BaseMapper finished: batches=" + std::to_string(stats.batches) +
           ", reads=" + std::to_string(stats.reads_processed) +
           ", mapped=" + std::to_string(stats.reads_mapped) +
           ", unmapped=" + std::to_string(stats.reads_unmapped));
  if (config_.result_writer) {
    LOG_INFO("Results written: " + std::to_string(stats.results_written) +
             " (primary=" + std::to_string(stats.primary_alignments) +
             ", secondary=" + std::to_string(stats.secondary_alignments) + ")");
  }
  return stats;
}

void BaseMapper::load_batch(BaseBatchBuffer& batch) {
  batch.num_reads = 0;
  const std::size_t capacity = batch.reads.size();
  while (batch.num_reads < capacity) {
    if (!provider_.next(batch.reads[batch.num_reads])) break;
    ++batch.num_reads;
  }
}

void BaseMapper::process_batch(BaseBatchBuffer& batch) {
  LOG_INFO("Processing batch: " + std::to_string(batch.num_reads) + " reads");

  executor_->parallel_for(0, batch.num_reads, 1,
                          [&](std::size_t i) { process_read(batch, i); });

  std::size_t total_hits = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    total_hits += batch.map_results[i].total_seed_hits;
  }
  LOG_INFO("Batch complete: total_hits=" + std::to_string(total_hits));
}

void BaseMapper::process_read(BaseBatchBuffer& batch, std::size_t index) {
  const auto t_start = std::chrono::high_resolution_clock::now();

  const auto& read = batch.reads[index];
  auto& seeds = batch.seeds[index];
  auto& hits = batch.seed_hits[index];

  ChainResult last_chains;
  std::size_t chunks_processed = 1;
  std::size_t total_hits = 0;

  if (config_.chunk_bp == 0 || read.sequence.size() <= config_.chunk_bp) {
    /* Whole-read path: extract once, look up, chain. */
    seeds = piru::base::extract_minimizers(read.sequence, config_.seeder);
    components_.lookup.lookup(seeds, hits);
    last_chains = components_.chainer->chain(hits);
    total_hits = hits.size();
  } else {
    /* Chunked path: extend prefix by chunk_bp each iteration; re-seed
     * cumulatively (bounded by max_chunks). After each chunk: chain,
     * check decision, early-exit if confident. */
    const std::size_t total_bases = read.sequence.size();
    const std::size_t max_chunks =
        config_.max_chunks > 0 ? config_.max_chunks : std::numeric_limits<std::size_t>::max();
    std::size_t chunk_idx = 0;
    while (chunk_idx < max_chunks) {
      const std::size_t chunk_end = std::min((chunk_idx + 1) * config_.chunk_bp, total_bases);
      std::string_view prefix(read.sequence.data(), chunk_end);

      seeds = piru::base::extract_minimizers(prefix, config_.seeder);
      hits.clear();
      components_.lookup.lookup(seeds, hits);

      ChainResult chunk_chains = components_.chainer->chain(hits);
      total_hits = hits.size();
      ++chunk_idx;

      bool decided = checkMappingDecision(
          chunk_chains.chains, config_.map_w_bestq, config_.map_w_bestmq, config_.map_w_bestmc,
          config_.map_w_threshold, config_.map_sc_min_anchors, config_.map_sc_ratio_lo,
          config_.map_sc_ratio_hi);

      last_chains = std::move(chunk_chains);

      if (!config_.no_early_exit && decided) break;
      if (chunk_end >= total_bases) break;
    }
    chunks_processed = chunk_idx;
  }

  hits.clear();
  hits.shrink_to_fit();

  /* Build map result from final chain. */
  ReadMapResult& result = batch.map_results[index];
  result.mappings.clear();
  result.total_seed_hits = total_hits;
  result.expanded_anchor_count = last_chains.expanded_anchor_count;
  result.chunks_processed = chunks_processed;
  const auto t_end = std::chrono::high_resolution_clock::now();
  result.processing_time_sec = std::chrono::duration<double>(t_end - t_start).count();

  for (const auto& chain : last_chains.chains) {
    Mapping mapping;
    mapping.anchors = chain.anchors;
    mapping.chain_score = chain.score;
    mapping.path_id = chain.path_id;
    mapping.coord_space = chain.coord_space;
    result.mappings.push_back(std::move(mapping));
  }

  /* Per-mapping mapq. */
  const int min_chain_sc = 15;
  for (std::size_t i = 0; i < result.mappings.size(); ++i) {
    double sec = (i + 1 < result.mappings.size()) ? result.mappings[i + 1].chain_score : 0.0;
    result.mappings[i].mapq =
        computeMapq(result.mappings[i].chain_score, sec, result.mappings[i].anchors.size(),
                    min_chain_sc);
  }

  /* Mapping decision applied to final chains. */
  float final_standout = 0.0f;
  DecisionPath decision_path = DecisionPath::kUnmapped;
  result.is_mapped = checkMappingDecision(
      last_chains.chains, config_.map_w_bestq, config_.map_w_bestmq, config_.map_w_bestmc,
      config_.map_w_threshold, config_.map_sc_min_anchors, config_.map_sc_ratio_lo,
      config_.map_sc_ratio_hi, &final_standout, &decision_path);
  result.standout = final_standout;

  /* EMA: only multi-chain standout-accepted chains contribute. Single-chain
   * gate uses its own ratio test, not the score distribution. */
  if (result.is_mapped && result.mappings.size() > 1 && decision_path == DecisionPath::kStandout) {
    recordAcceptedChain(result.mappings[0].chain_score, chunks_processed);
  }

  /* Fallback: multi-chain reads that failed standout. Threshold =
   * max(ema_at_ck, floor). Score-only, deterministic. */
  if (!result.is_mapped && result.mappings.size() > 1) {
    if (result.mappings[0].chain_score >= getFallbackThreshold(chunks_processed)) {
      result.is_mapped = true;
      decision_path = DecisionPath::kFallback;
    }
  }
  result.decision_path = decision_path;
}

BaseMapperStats BaseMapper::output_batch(const BaseBatchBuffer& batch) const {
  BaseMapperStats stats;

  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    const auto& read = batch.reads[i];
    const auto& map_result = batch.map_results[i];

    const bool is_mapped = map_result.mapped();
    if (is_mapped) {
      ++stats.reads_mapped;
    } else {
      ++stats.reads_unmapped;
    }

    double pri_score = is_mapped ? map_result.mappings[0].chain_score : 0.0;
    std::size_t pri_anchors = is_mapped ? map_result.mappings[0].anchors.size() : 0;
    LOG_DEBUG("DECISION\t" + read.id + "\t" + (is_mapped ? "MAP" : "UNMAP") + "\t" +
              "score=" + std::to_string(static_cast<int>(pri_score)) + "\t" +
              "anchors=" + std::to_string(pri_anchors));

    piru::io::ResultWriter* writer = config_.result_writer;
    if (!config_.per_file_writers.empty()) {
      auto it = config_.per_file_writers.find(read.source_file);
      if (it != config_.per_file_writers.end()) writer = it->second;
    }
    if (writer) {
      // base mode has no raw signal length; pass sequence length so
      // downstream eval has SOMETHING for read-length tagging.
      writer->write(map_result, read.id, read.sequence.size());
      if (is_mapped) {
        ++stats.primary_alignments;
        std::size_t sec_count = map_result.mappings.size() > 1 ? map_result.mappings.size() - 1 : 0;
        stats.secondary_alignments += sec_count;
        stats.results_written += 1 + sec_count;
      }
    }
  }

  return stats;
}

}  // namespace piru::base::mapping
