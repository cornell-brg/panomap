/**
 * batch_mapper.cpp
 *
 * Orchestrates read mapping: signal processing, seed lookup, chaining,
 * mapping decision, and result output. Processes reads in batches.
 *
 * Related:
 *  - batch_mapper.hpp
 *  - path_chainer.cpp, graph_chainer.cpp, sort_chainer.cpp (chaining backends)
 *  - gaf_writer.cpp (result output)
 *
 * SPDX-License-Identifier: MIT
 */

#include "signal/mapping/batch_mapper.hpp"

#include "signal/diff_filter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <utility>

#include "core/mapping/pan_chainer.hpp"
#include "core/mapping/path_chainer.hpp"
#include "core/mapping/sort_chainer.hpp"
#include "core/util/logging.hpp"
#include "core/util/trace.hpp"

namespace panomap::mapping {

namespace {

// Compute mapq for a chain.
// (1) Absolute score strength (0 to score_range)
// (2) Relative standout vs next chain (0 to standout_range)
// (3) Anchor count robustness (caps everything)
// standout_ratio controls the split: standout gets ratio*60, score gets (1-ratio)*60.
// RH2-style mapq: log-based, with anchor count and secondary suppression.
// Matches mm_set_mapq() in hit.c (non-DTW mode).
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

// Compute event_span (query) and ref_span from a chain's anchors.
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

// Unified mapping decision: is the best chain confident enough?
// Used for both early exit (per-chunk) and final map/unmap decision.
//
// Single-chain reads: must have anchors >= sc_min_anchors AND event/ref ratio in [sc_ratio_lo, sc_ratio_hi].
// Multi-chain reads: weighted standout of best chain vs others must exceed w_threshold.
//
// out_standout: if non-null, receives the weighted standout value.
// out_path: if non-null, receives the decision path (kGate/kStandout/kUnmapped). Fallback
//           decisions are applied later in the caller.
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

  // Always compute standout so we can report it (even for single-chain reads).
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

  // Single-chain gate: anchor count + event/ref ratio
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

  // Multi-chain: weighted standout threshold
  if (weighted >= w_threshold) {
    if (out_path) *out_path = DecisionPath::kStandout;
    return true;
  }
  return false;
}

}  // namespace

void SeedLookup::lookup(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& out_hits) const {
  if (!store_) return;
  out_hits.clear();
  out_hits.reserve(seeds.seeds.size());
  for (const auto& seed : seeds.seeds) {
    auto hits = store_->lookup(seed.hash);
    if (hits.count == 0) continue;
    if (hits.count > freq_threshold_) continue;  // skip overly frequent seeds
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
    // Per-read hit cap: stop collecting if we already have enough hits.
    // Prevents OOM when many seeds match high-frequency index entries.
    if (max_total_hits_ > 0 && out_hits.size() >= max_total_hits_) break;
  }
}

void BatchBuffer::resize(std::size_t capacity) {
  raw_reads.resize(capacity);
  normalized.resize(capacity);
  tokenized.resize(capacity);
  seeds.resize(capacity);
  seed_hits.resize(capacity);
  map_results.resize(capacity);
  num_reads = 0;
}

void BatchBuffer::clear() {
  for (std::size_t i = 0; i < num_reads; ++i) {
    raw_reads[i] = io::RawRead{};
    normalized[i] = signal::NormalizedSignal{};
    tokenized[i].tokens.clear();
    seeds[i].seeds.clear();
    seed_hits[i].clear();
    map_results[i] = ReadMapResult{};
  }
  num_reads = 0;
}

BatchMapper::BatchMapper(io::ReadProvider& provider, BatchMapperConfig config, std::ostream& output)
    : config_(std::move(config)),
      provider_(provider),
      executor_(concurrency::make_executor(config_.num_threads)),
      components_(create_components()),
      output_(output) {
  // Pre-size per-chunk-depth EMA vector. Index 0 unused, 1..max_chunks active.
  // EMA starts at the floor value so threshold = max(ema, floor) begins at floor.
  const std::size_t mc = config_.event_pipeline_config.max_chunks;
  const std::size_t sz = (mc > 0 ? mc : 20) + 1;
  ema_score_per_ck_.assign(sz, config_.map_fallback_floor_score);
}

void BatchMapper::recordAcceptedChain(double score, std::size_t chunks_processed) {
  if (!config_.map_fallback_adaptive) return;
  if (chunks_processed == 0 || chunks_processed >= ema_score_per_ck_.size()) return;
  std::lock_guard<std::mutex> lock(adaptive_mutex_);
  const double a = config_.map_fallback_alpha;
  ema_score_per_ck_[chunks_processed] =
      a * score + (1.0 - a) * ema_score_per_ck_[chunks_processed];
}

double BatchMapper::getFallbackThreshold(std::size_t chunks_processed) const {
  const double floor = config_.map_fallback_floor_score;
  if (!config_.map_fallback_adaptive) return floor;
  if (chunks_processed == 0 || chunks_processed >= ema_score_per_ck_.size()) return floor;
  std::lock_guard<std::mutex> lock(adaptive_mutex_);
  // Conservative: max of learned mean and absolute floor. Floor prevents
  // EMA drift from dragging threshold into noise territory.
  return std::max(ema_score_per_ck_[chunks_processed], floor);
}

PipelineComponents BatchMapper::create_components() const {
  PipelineComponents comps;

  /* Create unified event pipeline (event detection + normalization) */
  comps.event_pipeline = signal::make_event_pipeline(config_.event_pipeline_config);

  comps.tokenizer = signal::make_tokenizer(config_.tokenizer_config);
  comps.seed_extractor = signal::make_seed_extractor(config_.seed_config);
  comps.seed_store = config_.seed_store;
  comps.graph_store = config_.graph_store;
  if (!comps.seed_store) {
    throw std::runtime_error("BatchMapper requires a SeedStore for lookup");
  }

  /* Create chainer backend */
  if (config_.chainer_backend == "path-chain") {
    if (!config_.linearization_coords) {
      throw std::runtime_error("PathChainer requires linearization_coords");
    }
    if (!config_.path_lengths) {
      throw std::runtime_error("PathChainer requires path_lengths for bounds checking");
    }
    auto path_config = PathChainerConfig::from_parsed(config_.chainer_parsed);
    path_config.merge_anchors = config_.enable_anchor_merge;
    path_config.pore_k = config_.pore_k;
    // Normalize gap penalty by anchor span (matching RH2: scale * 0.01 * (e + k - 1))
    if (config_.pore_k > 0 && config_.seed_config.k > 0) {
      float span = static_cast<float>(config_.seed_config.k + config_.pore_k - 1);
      path_config.chn_pen_gap = path_config.chn_pen_gap * 0.01f * span;
      path_config.chn_pen_skip = path_config.chn_pen_skip * 0.01f * span;
    }
    comps.chainer = std::make_unique<PathChainer>(path_config, *config_.linearization_coords,
                                                  *config_.path_lengths, config_.node_1d_coords,
                                                  config_.component_ids);
  } else if (config_.chainer_backend == "sort-chain") {
    if (!config_.node_1d_coords) {
      throw std::runtime_error(
          "SortChainer requires node_1d_coords (use --compute-1d-sort or --1d-coords-file at index "
          "time)");
    }
    auto sort_config = SortChainerConfig::from_parsed(config_.chainer_parsed);
    sort_config.pore_k = config_.pore_k;
    if (config_.pore_k > 0 && config_.seed_config.k > 0) {
      float span = static_cast<float>(config_.seed_config.k + config_.pore_k - 1);
      sort_config.chn_pen_gap = sort_config.chn_pen_gap * 0.01f * span;
      sort_config.chn_pen_skip = sort_config.chn_pen_skip * 0.01f * span;
    }
    std::vector<std::uint32_t> bp_lens(config_.graph_store->nodeCount());
    for (std::size_t i = 0; i < bp_lens.size(); ++i)
      bp_lens[i] = static_cast<std::uint32_t>(config_.graph_store->sequenceLen(i));
    static const std::vector<std::uint32_t> kEmptyCC;
    const auto& cc_ids = config_.component_ids ? *config_.component_ids : kEmptyCC;
    comps.chainer = std::make_unique<SortChainer>(sort_config, *config_.node_1d_coords,
                                                  std::move(bp_lens), cc_ids);
  } else if (config_.chainer_backend == "pan-chain") {
    if (!config_.node_1d_coords) {
      throw std::runtime_error(
          "PanChainer requires node_1d_coords (use --compute-1d-sort or --1d-coords-file)");
    }
    if (!config_.linearization_coords) {
      throw std::runtime_error("PanChainer requires linearization_coords");
    }
    if (!config_.path_lengths) {
      throw std::runtime_error("PanChainer requires path_lengths");
    }
    auto pan_config = PanChainerConfig::from_parsed(config_.chainer_parsed);
    pan_config.pore_k = config_.pore_k;
    if (config_.pore_k > 0 && config_.seed_config.k > 0) {
      float span = static_cast<float>(config_.seed_config.k + config_.pore_k - 1);
      pan_config.chn_pen_gap = pan_config.chn_pen_gap * 0.01f * span;
      pan_config.chn_pen_skip = pan_config.chn_pen_skip * 0.01f * span;
    }
    comps.chainer = std::make_unique<PanChainer>(
        pan_config, *config_.node_1d_coords, *config_.linearization_coords, *config_.path_lengths);
  } else {
    throw std::runtime_error("Unknown chainer backend: " + config_.chainer_backend);
  }
  const std::size_t freq_threshold = comps.seed_store->frequency_threshold();
  comps.lookup = SeedLookup(comps.seed_store, freq_threshold, config_.max_total_hits);

  /* Log pipeline configuration */
  LOG_DEBUG("Pipeline: " + comps.chainer->name() + " chaining");

  return comps;
}

BatchMapperStats BatchMapper::process_all() {
  BatchMapperStats stats;
  BatchBuffer batch;
  batch.resize(config_.batch_capacity_reads);

  LOG_DEBUG(
      "BatchMapper starting: batch_capacity_reads=" + std::to_string(config_.batch_capacity_reads) +
      ", threads=" + std::to_string(executor_->max_concurrency()) + ", tokenizer=" +
      components_.tokenizer->name() + ", seeds=" + components_.seed_extractor->name());

  while (true) {
    load_batch(batch);
    if (batch.num_reads == 0) break;

    ++stats.batches;
    stats.reads_processed += batch.num_reads;
    process_batch(batch);

    /* Accumulate output stats */
    auto batch_stats = output_batch(batch);
    stats.reads_mapped += batch_stats.reads_mapped;
    stats.reads_unmapped += batch_stats.reads_unmapped;
    stats.results_written += batch_stats.results_written;
    stats.primary_alignments += batch_stats.primary_alignments;
    stats.secondary_alignments += batch_stats.secondary_alignments;

    batch.clear();
  }

  /* Log detailed summary */
  LOG_INFO("BatchMapper finished: batches=" + std::to_string(stats.batches) +
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

void BatchMapper::load_batch(BatchBuffer& batch) {
  batch.num_reads = 0;
  const std::size_t capacity = batch.raw_reads.size();
  io::RawRead read;
  while (batch.num_reads < capacity && provider_.get_next(read)) {
    batch.raw_reads[batch.num_reads] = std::move(read);
    ++batch.num_reads;
  }
}

void BatchMapper::process_batch(BatchBuffer& batch) {
  /* Calculate estimated memory before processing */
  std::size_t est_mem_mb = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    est_mem_mb += batch.raw_reads[i].len_raw_signal * sizeof(int16_t) / 1024 / 1024;
  }
  LOG_INFO("Processing batch: " + std::to_string(batch.num_reads) + " reads, ~" +
           std::to_string(est_mem_mb) + " MB raw signal");

  executor_->parallel_for(0, batch.num_reads, 1, [&](std::size_t i) { process_read(batch, i); });

  /* Summarize seed hits from per-read results */
  std::size_t total_hits = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    total_hits += batch.map_results[i].total_seed_hits;
  }
  LOG_INFO("Batch complete: total_hits=" + std::to_string(total_hits));
}

void BatchMapper::process_read(BatchBuffer& batch, std::size_t index) {
  const auto t_start = std::chrono::high_resolution_clock::now();

  const auto& read = batch.raw_reads[index];
  auto& all_hits = batch.seed_hits[index];
  all_hits.clear();

  const std::size_t chunk_size = config_.event_pipeline_config.chunk_size;
  std::size_t chunks_processed = 1;
  ChainResult last_chains;  // holds the final chain result (from last chunk or single chain)

  if (chunk_size == 0) {
    /* Non-chunked path: process entire signal at once (original behavior) */
    batch.normalized[index] = components_.event_pipeline->process(read);
    signal::apply_diff_filter(batch.normalized[index], config_.diff_filter);
    batch.tokenized[index] = components_.tokenizer->quantize(batch.normalized[index]);

    batch.seeds[index] = components_.seed_extractor->extract(batch.tokenized[index]);
    components_.lookup.lookup(batch.seeds[index], all_hits);
  } else {
    /* Chunked path: process signal in chunks, accumulate seed hits.
     * After each chunk: chain, filter survivors, check early exit.
     * Like RH2: independent DSP per chunk, cumulative hits for chaining,
     * early exit when confident mapping found. */

    /* Convert raw ADC to picoamps (once, full signal) */
    const float raw_unit = (read.digitisation == 0.0f) ? 1.0f : (read.range / read.digitisation);
    std::vector<float> pA;
    pA.reserve(read.raw_signal.size());
    for (auto value : read.raw_signal) {
      float pa = (static_cast<float>(value) + read.offset) * raw_unit;
      if (pa > 30.0f && pa < 200.0f) {
        pA.push_back(pa);
      }
    }

    const std::size_t total_samples = pA.size();
    const std::size_t max_chunks = config_.event_pipeline_config.max_chunks;
    signal::NormState norm_state;
    std::size_t chunk_idx = 0;
    std::uint32_t query_offset = 0;  // cumulative event count for query position offset

    for (std::size_t chunk_start = 0; chunk_start < total_samples; chunk_start += chunk_size) {
      // Max chunks limit (0 = unlimited)
      if (max_chunks > 0 && chunk_idx >= max_chunks) break;

      std::size_t chunk_len = std::min(chunk_size, total_samples - chunk_start);

      /* Per-chunk DSP: event detection, quantization, seed extraction */
      auto events =
          components_.event_pipeline->process_chunk(pA.data() + chunk_start, chunk_len, norm_state);
      signal::apply_diff_filter(events, config_.diff_filter);
      auto tokens = components_.tokenizer->quantize(events);

      auto seeds = components_.seed_extractor->extract(tokens);

      /* Lookup and accumulate hits, offsetting query positions by cumulative event count */
      std::vector<NodeAnchor> chunk_hits;
      components_.lookup.lookup(seeds, chunk_hits);
      for (auto& hit : chunk_hits) {
        hit.read_pos += query_offset;
      }
      all_hits.insert(all_hits.end(), chunk_hits.begin(), chunk_hits.end());

      // like RH2's reg->offset
      query_offset += static_cast<std::uint32_t>(events.samples.size());

      /* Trace: per-chunk pipeline dumps */
      const auto& rid = read.read_id;
      PANOMAP_TRACE_DUMP(trace::kSignal, rid, {
        std::ofstream ofs(trace::trace_path("1_pA", rid, chunk_idx));
        for (std::size_t i = 0; i < chunk_len; ++i) ofs << pA[chunk_start + i] << "\n";
      });
      PANOMAP_TRACE_DUMP(trace::kNorm, rid, {
        std::ofstream ofs(trace::trace_path("1b_norm", rid, chunk_idx));
        for (const auto& v : events.debug_normalized) ofs << v << "\n";
      });
      PANOMAP_TRACE_DUMP(trace::kEvents, rid, {
        std::ofstream ofs(trace::trace_path("2_events", rid, chunk_idx));
        for (const auto& v : events.samples) ofs << v << "\n";
      });
      PANOMAP_TRACE_DUMP(trace::kTokens, rid, {
        std::ofstream ofs(trace::trace_path("3_tokens", rid, chunk_idx));
        for (const auto& v : tokens.tokens) ofs << v << "\n";
      });
      PANOMAP_TRACE_DUMP(trace::kSeeds, rid, {
        std::ofstream ofs(trace::trace_path("4_seeds", rid, chunk_idx));
        for (const auto& s : seeds.seeds) ofs << s.hash << "\t" << s.position << "\n";
      });
      PANOMAP_TRACE_DUMP(trace::kHits, rid, {
        std::ofstream ofs(trace::trace_path("5_hits", rid, chunk_idx));
        ofs << "chunk_new=" << chunk_hits.size() << "\n";
        for (const auto& h : chunk_hits)
          ofs << h.node_id << "\t" << h.offset << "\t" << h.read_pos << "\t" << h.span << "\n";
      });

      /* Chain all accumulated hits, keep only survivors */
      ChainResult chunk_chains = components_.chainer->chain(all_hits);

      PANOMAP_TRACE_DUMP(trace::kChains, rid, {
        std::ofstream ofs(trace::trace_path("6_anchors", rid, chunk_idx));
        ofs << "# Anchor cloud: query_pos\t1d_ref_coord\tnode_id\toffset\tspan\n";
        for (const auto& h : all_hits) {
          double ref_1d = 0.0;
          if (config_.node_1d_coords && h.node_id < config_.node_1d_coords->size())
            ref_1d = (*config_.node_1d_coords)[h.node_id] + static_cast<double>(h.offset);
          ofs << h.read_pos << "\t" << ref_1d << "\t" << h.node_id << "\t" << h.offset << "\t"
              << h.span << "\n";
        }
        ofs << "# Chains: " << chunk_chains.chains.size() << "\n";
        for (std::size_t ci = 0; ci < chunk_chains.chains.size(); ++ci) {
          const auto& chain = chunk_chains.chains[ci];
          ofs << "# chain " << ci << " score=" << chain.score << " anchors=" << chain.anchors.size()
              << "\n";
          for (const auto& a : chain.anchors) {
            ofs << a.read_pos << "\t" << a.ref_coord << "\t" << a.node_id << "\t" << a.offset
                << "\t" << a.length << "\tCHAIN=" << ci << "\n";
          }
        }
      });

      if (!chunk_chains.used_inputs.empty()) {
        std::vector<NodeAnchor> survivors;
        survivors.reserve(all_hits.size());
        for (std::size_t i = 0; i < all_hits.size(); ++i) {
          if (i < chunk_chains.used_inputs.size() && chunk_chains.used_inputs[i]) {
            survivors.push_back(all_hits[i]);
          }
        }
        LOG_DEBUG("SURVIVORS\t" + read.read_id + "\tchunk=" + std::to_string(chunk_idx) +
                  "\tbefore=" + std::to_string(all_hits.size()) +
                  "\tafter=" + std::to_string(survivors.size()) +
                  "\tchains=" + std::to_string(chunk_chains.chains.size()));
        all_hits = std::move(survivors);
      }

      /* Check mapping decision (early exit or last chunk) */
      bool decided = checkMappingDecision(chunk_chains.chains, config_.map_w_bestq,
                                          config_.map_w_bestmq, config_.map_w_bestmc,
                                          config_.map_w_threshold, config_.map_sc_min_anchors,
                                          config_.map_sc_ratio_lo, config_.map_sc_ratio_hi);

      if (!config_.no_early_exit && decided) {
        last_chains = std::move(chunk_chains);
        ++chunk_idx;
        break;
      }

      last_chains = std::move(chunk_chains);
      ++chunk_idx;
    }
    chunks_processed = chunk_idx;
  }

  /* Use chunk chains directly (no re-chain). For non-chunked path, chain once. */
  ChainResult& chain_result =
      last_chains.chains.empty()
          ? (last_chains = components_.chainer->chain(all_hits))  // non-chunked: first chain
          : last_chains;                                          // chunked: last chunk's chains

  const std::size_t total_hits = all_hits.size();
  all_hits.clear();
  all_hits.shrink_to_fit();

  /* Build map result from chain result */
  ReadMapResult& result = batch.map_results[index];
  result.mappings.clear();
  result.total_seed_hits = total_hits;
  result.expanded_anchor_count = chain_result.expanded_anchor_count;
  result.chunks_processed = chunks_processed;
  const auto t_end = std::chrono::high_resolution_clock::now();
  result.processing_time_sec = std::chrono::duration<double>(t_end - t_start).count();

  for (const auto& chain : chain_result.chains) {
    Mapping mapping;
    mapping.anchors = chain.anchors;
    mapping.chain_score = chain.score;
    mapping.path_id = chain.path_id;
    mapping.coord_space = chain.coord_space;
    result.mappings.push_back(std::move(mapping));
  }

  /* Compute mapq for all mappings */
  const int min_chain_sc = 15;
  for (std::size_t i = 0; i < result.mappings.size(); ++i) {
    double sec = (i + 1 < result.mappings.size()) ? result.mappings[i + 1].chain_score : 0.0;
    result.mappings[i].mapq =
        computeMapq(result.mappings[i].chain_score, sec, result.mappings[i].anchors.size(),
                    min_chain_sc);
  }

  /* Mapping decision: same function used for early exit, applied to final chains */
  float final_standout = 0.0f;
  DecisionPath decision_path = DecisionPath::kUnmapped;
  result.is_mapped = checkMappingDecision(
      chain_result.chains, config_.map_w_bestq, config_.map_w_bestmq, config_.map_w_bestmc,
      config_.map_w_threshold, config_.map_sc_min_anchors, config_.map_sc_ratio_lo,
      config_.map_sc_ratio_hi, &final_standout, &decision_path);
  result.standout = final_standout;

  /* Record accepted chain score at this chunk depth.
   * Only records multi-chain accepts; single-chain uses its own gate, not EMA. */
  if (result.is_mapped && !result.mappings.empty() && result.mappings.size() > 1) {
    recordAcceptedChain(result.mappings[0].chain_score, chunks_processed);
  }

  /* Fallback: multi-chain reads that failed standout. Score-only check.
   * Threshold = max(ema_mean, floor). Floor prevents EMA drift to noise territory. */
  if (!result.is_mapped && result.mappings.size() > 1) {
    double thresh = getFallbackThreshold(chunks_processed);
    const auto& best = result.mappings[0];
    if (best.chain_score >= thresh) {
      result.is_mapped = true;
      decision_path = DecisionPath::kFallback;
    }
  }
  result.decision_path = decision_path;

}

BatchMapperStats BatchMapper::output_batch(const BatchBuffer& batch) const {
  BatchMapperStats stats;

  /* Try to get path names from graph store */
  const index::AdjListGraphStore* adj_store =
      dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);

  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    const auto& read = batch.raw_reads[i];
    const auto& map_result = batch.map_results[i];

    /* Track mapped/unmapped */
    const bool is_mapped = map_result.mapped();
    if (is_mapped) {
      ++stats.reads_mapped;
    } else {
      ++stats.reads_unmapped;
    }

    /* Per-read decision log (to stderr when verbose) */
    double pri_score = is_mapped ? map_result.mappings[0].chain_score : 0.0;
    std::size_t pri_anchors = is_mapped ? map_result.mappings[0].anchors.size() : 0;
    LOG_DEBUG("DECISION\t" + read.read_id + "\t" + (is_mapped ? "MAP" : "UNMAP") + "\t" +
              "chunks=" + std::to_string(map_result.chunks_processed) + "\t" +
              "score=" + std::to_string(static_cast<int>(pri_score)) + "\t" +
              "anchors=" + std::to_string(pri_anchors));

    /* Write to result file if configured.
     * Per-file routing: look up writer by source_file, fall back to single writer. */
    io::ResultWriter* writer = config_.result_writer;
    if (!config_.per_file_writers.empty()) {
      auto it = config_.per_file_writers.find(read.source_file);
      if (it != config_.per_file_writers.end()) {
        writer = it->second;
      }
    }
    if (writer) {
      writer->write(map_result, read.read_id, read.len_raw_signal);
      if (is_mapped) {
        ++stats.primary_alignments;
        std::size_t sec_count = map_result.mappings.size() > 1 ? map_result.mappings.size() - 1 : 0;
        stats.secondary_alignments += sec_count;
        stats.results_written += 1 + sec_count;
      }
    }

    // Suppress unused variable warnings
    (void)adj_store;
  }

  return stats;
}

}  // namespace panomap::mapping
