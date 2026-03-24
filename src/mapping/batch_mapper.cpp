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

#include "mapping/batch_mapper.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <utility>

#include "mapping/graph_chainer.hpp"
#include "mapping/path_chainer.hpp"
#include "mapping/sort_chainer.hpp"
#include "util/logging.hpp"
#include "util/trace.hpp"

namespace piru::mapping {

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
  fuzzy_quantized.resize(capacity);
  seeds.resize(capacity);
  seed_hits.resize(capacity);
  map_results.resize(capacity);
  num_reads = 0;
}

void BatchBuffer::clear() {
  for (std::size_t i = 0; i < num_reads; ++i) {
    raw_reads[i] = io::RawRead{};
    normalized[i] = signal::NormalizedSignal{};
    fuzzy_quantized[i].tokens.clear();
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
      output_(output) {}

PipelineComponents BatchMapper::create_components() const {
  PipelineComponents comps;

  /* Create unified event pipeline (event detection + normalization) */
  comps.event_pipeline = signal::make_event_pipeline(config_.event_pipeline_config);

  comps.fuzzy_quantizer = signal::make_fuzzy_quantizer(config_.fuzzy_config);
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
    comps.chainer = std::make_unique<PathChainer>(path_config, *config_.linearization_coords,
                                                  *config_.path_lengths);
  } else if (config_.chainer_backend == "graph-chain") {
    if (!config_.linearization_coords) {
      throw std::runtime_error("GraphChainer requires linearization_coords");
    }
    if (!config_.path_lengths) {
      throw std::runtime_error("GraphChainer requires path_lengths for bounds checking");
    }
    comps.chainer = std::make_unique<GraphChainer>(*config_.linearization_coords,
                                                   *config_.path_lengths);
  } else if (config_.chainer_backend == "sort-chain") {
    if (!config_.node_1d_coords) {
      throw std::runtime_error("SortChainer requires node_1d_coords (use --compute-1d-sort or --1d-coords-file at index time)");
    }
    SortChainerConfig sort_config;
    sort_config.pore_k = config_.pore_k;
    comps.chainer = std::make_unique<SortChainer>(sort_config, *config_.node_1d_coords);
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
      ", threads=" + std::to_string(executor_->max_concurrency()) + ", fuzzy=" +
      components_.fuzzy_quantizer->name() + ", seeds=" + components_.seed_extractor->name());

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

  /* Calculate seed hits memory */
  std::size_t total_hits = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    total_hits += batch.seed_hits[i].size();
  }
  const std::size_t hits_mem_mb = total_hits * sizeof(NodeAnchor) / 1024 / 1024;
  LOG_INFO("Batch complete: total_hits=" + std::to_string(total_hits) + " (~" +
           std::to_string(hits_mem_mb) + " MB)");
}

void BatchMapper::process_read(BatchBuffer& batch, std::size_t index) const {
  const auto t_start = std::chrono::high_resolution_clock::now();

  const auto& read = batch.raw_reads[index];
  auto& all_hits = batch.seed_hits[index];
  all_hits.clear();

  const std::size_t chunk_size = config_.event_pipeline_config.chunk_size;
  std::size_t chunks_processed = 1;  // for mapping decision score scaling

  if (chunk_size == 0) {
    /* Non-chunked path: process entire signal at once (original behavior) */
    batch.normalized[index] = components_.event_pipeline->process(read);
    batch.fuzzy_quantized[index] = components_.fuzzy_quantizer->quantize(batch.normalized[index]);
    batch.seeds[index] = components_.seed_extractor->extract(batch.fuzzy_quantized[index]);
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
      pA.push_back((static_cast<float>(value) + read.offset) * raw_unit);
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
      auto events = components_.event_pipeline->process_chunk(
          pA.data() + chunk_start, chunk_len, norm_state);
      auto tokens = components_.fuzzy_quantizer->quantize(events);
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
      PIRU_TRACE_DUMP(trace::kSignal, rid, {
        std::ofstream ofs(trace::trace_path("1_pA", rid, chunk_idx));
        for (std::size_t i = 0; i < chunk_len; ++i) ofs << pA[chunk_start + i] << "\n";
      });
      PIRU_TRACE_DUMP(trace::kNorm, rid, {
        std::ofstream ofs(trace::trace_path("1b_norm", rid, chunk_idx));
        for (const auto& v : events.debug_normalized) ofs << v << "\n";
      });
      PIRU_TRACE_DUMP(trace::kEvents, rid, {
        std::ofstream ofs(trace::trace_path("2_events", rid, chunk_idx));
        for (const auto& v : events.samples) ofs << v << "\n";
      });
      PIRU_TRACE_DUMP(trace::kTokens, rid, {
        std::ofstream ofs(trace::trace_path("3_tokens", rid, chunk_idx));
        for (const auto& v : tokens.tokens) ofs << v << "\n";
      });
      PIRU_TRACE_DUMP(trace::kSeeds, rid, {
        std::ofstream ofs(trace::trace_path("4_seeds", rid, chunk_idx));
        for (const auto& s : seeds.seeds) ofs << s.hash << "\t" << s.position << "\n";
      });
      PIRU_TRACE_DUMP(trace::kHits, rid, {
        std::ofstream ofs(trace::trace_path("5_hits", rid, chunk_idx));
        ofs << "chunk_new=" << chunk_hits.size() << "\n";
        for (const auto& h : chunk_hits)
          ofs << h.node_id << "\t" << h.offset << "\t" << h.read_pos << "\t" << h.span << "\n";
      });

      /* Chain all accumulated hits, keep only survivors */
      ChainResult chunk_chains = components_.chainer->chain(all_hits);
      if (!chunk_chains.used_inputs.empty()) {
        std::vector<NodeAnchor> survivors;
        survivors.reserve(all_hits.size());
        for (std::size_t i = 0; i < all_hits.size(); ++i) {
          if (i < chunk_chains.used_inputs.size() && chunk_chains.used_inputs[i]) {
            survivors.push_back(all_hits[i]);
          }
        }
        all_hits = std::move(survivors);
      }

      /* Early exit: if best chain has high enough MAPQ, stop processing chunks */
      if (config_.map_threshold > 0.0f && !chunk_chains.chains.empty()) {
        double best_cs = chunk_chains.chains[0].score;
        double sec_cs = (chunk_chains.chains.size() > 1) ? chunk_chains.chains[1].score : 0.0;
        double mapq = (best_cs <= 0.0) ? 0.0
            : (sec_cs <= 0.0) ? 60.0
            : std::clamp(40.0 * (1.0 - sec_cs / best_cs), 0.0, 60.0);
        // RH2 early exit: single chain with mapq >= min_mapq (default 2)
        if (chunk_chains.chains.size() == 1 && mapq >= 2.0) break;
      }

      ++chunk_idx;
    }
    chunks_processed = chunk_idx;
  }

  /* ROI anchor filtering */
  if (config_.roi_filter_anchors && config_.roi_nodes) {
    const auto& roi = *config_.roi_nodes;
    all_hits.erase(std::remove_if(all_hits.begin(), all_hits.end(),
                                   [&roi](const NodeAnchor& a) {
                                     return roi.count(a.node_id) == 0;
                                   }),
                   all_hits.end());
  }

  /* Final chain with all accumulated hits */
  ChainResult chain_result = components_.chainer->chain(all_hits);

  /* Free seed hits */
  all_hits.clear();
  all_hits.shrink_to_fit();

  /* Build unified map result from chain result */
  ReadMapResult& result = batch.map_results[index];
  result.mappings.clear();
  result.expanded_anchor_count = chain_result.expanded_anchor_count;
  result.chunks_processed = chunks_processed;
  const auto t_end = std::chrono::high_resolution_clock::now();
  result.processing_time_sec =
      std::chrono::duration<double>(t_end - t_start).count();

  for (const auto& chain : chain_result.chains) {
    Mapping mapping;
    mapping.anchors = chain.anchors;
    mapping.chain_score = chain.score;
    mapping.path_id = chain.path_id;
    mapping.coord_space = chain.coord_space;
    result.mappings.push_back(std::move(mapping));
  }

  /* Mapping decision: RH2-style weighted score filter.
   * Computes how much the best chain stands out from the average.
   * If weighted score < threshold, declare unmapped. */
  if (config_.map_threshold > 0.0f && result.mappings.size() >= 1) {
    double best_score = result.mappings[0].chain_score;

    /* Compute MAPQ from best/secondary ratio */
    double secondary_score = (result.mappings.size() > 1)
        ? result.mappings[1].chain_score : 0.0;
    double best_mapq = (best_score <= 0.0) ? 0.0
        : (secondary_score <= 0.0) ? 60.0
        : std::clamp(40.0 * (1.0 - secondary_score / best_score), 0.0, 60.0);

    /* Per-chain mapq: each chain's mapq is based on its score vs the next chain.
     * Mean score and mean mapq across all chains. */
    double mean_score = 0.0, mean_mapq = 0.0;
    const auto& mappings = result.mappings;
    for (std::size_t mi = 0; mi < mappings.size(); ++mi) {
      mean_score += mappings[mi].chain_score;
      double sec = (mi + 1 < mappings.size()) ? mappings[mi + 1].chain_score : 0.0;
      double mq = (mappings[mi].chain_score <= 0.0) ? 0.0
          : (sec <= 0.0) ? 60.0
          : std::clamp(40.0 * (1.0 - sec / mappings[mi].chain_score), 0.0, 60.0);
      mean_mapq += mq;
    }
    mean_score /= static_cast<double>(mappings.size());
    mean_mapq /= static_cast<double>(mappings.size());

    /* Weighted components */
    float scaled_score_scale = config_.map_score_scale * std::sqrt(static_cast<float>(chunks_processed));
    float chunk_factor = std::max(0.5f, 1.0f - 1.0f / static_cast<float>(chunks_processed));
    float effective_w_abs = config_.map_w_abs * chunk_factor;
    float r_abs = std::min(static_cast<float>(best_score / scaled_score_scale), 1.0f);
    float r_bestq = (best_mapq > 0.0) ? std::min(static_cast<float>(best_mapq / 30.0), 1.0f) : 0.0f;
    // Clamp to [0,1] matching RH2
    float r_bestmq = (best_mapq > 0.0)
        ? std::clamp(static_cast<float>(1.0 - mean_mapq / best_mapq), 0.0f, 1.0f)
        : 0.0f;
    float r_bestmc = (best_score > 0.0)
        ? std::clamp(static_cast<float>(1.0 - mean_score / best_score), 0.0f, 1.0f)
        : 0.0f;

    float weighted = effective_w_abs * r_abs
                   + config_.map_w_bestq * r_bestq
                   + config_.map_w_bestmq * r_bestmq
                   + config_.map_w_bestmc * r_bestmc;

    PIRU_TRACE_DUMP(trace::kChains, read.read_id, {
      std::ofstream ofs(trace::trace_path("7_decision", read.read_id));
      ofs << "bs=" << best_score << "\tss=" << secondary_score
          << "\tms=" << mean_score << "\tbq=" << best_mapq
          << "\tmq=" << mean_mapq << "\tw=" << weighted
          << "\tn=" << result.mappings.size() << "\n";
    });

    if (weighted < config_.map_threshold) {
      result.mappings.clear();  // unmapped
    }
  }

  /* Compute per-mapping MAPQ */
  if (!result.mappings.empty()) {
    for (std::size_t i = 0; i < result.mappings.size(); ++i) {
      double sec = (i + 1 < result.mappings.size()) ? result.mappings[i + 1].chain_score : 0.0;
      double pri = result.mappings[i].chain_score;
      if (pri <= 0.0) {
        result.mappings[i].mapq = 0;
      } else if (sec <= 0.0) {
        result.mappings[i].mapq = 60;
      } else {
        result.mappings[i].mapq = static_cast<int>(
            std::clamp(40.0 * (1.0 - sec / pri), 0.0, 60.0) + 0.499);
      }
    }

    /* Secondary suppression: if many chains score similarly and primary
     * mapq is low, only keep the primary. */
    if (result.mappings.size() > 1 && result.mappings[0].mapq < 5) {
      double primary_score = result.mappings[0].chain_score;
      std::size_t n_similar = 0;
      for (const auto& m : result.mappings) {
        if (m.chain_score >= primary_score * 0.8) ++n_similar;
      }
      if (n_similar > 3) {
        result.mappings.resize(1);
      }
    }
  }

  /* ROI classification (if --roi is active) */
  if (config_.roi_nodes && !config_.classify_mode.empty() && result.mapped()) {
    bool above_threshold;

    if (config_.roi_filter_anchors) {
      /* --chain-target: classify by chain score threshold.
       * All anchors are ROI nodes, so overlap is meaningless. */
      double chain_score = result.mappings[0].chain_score;
      result.roi_overlap = chain_score;  // repurpose field for output
      above_threshold = chain_score >= config_.roi_score_threshold;
    } else {
      /* --chain-genome: classify by ROI node overlap fraction */
      const auto& anchors = result.mappings[0].anchors;
      const auto& roi = *config_.roi_nodes;

      double roi_bases = 0.0;
      double total_bases = 0.0;

      for (std::size_t i = 0; i + 1 < anchors.size(); ++i) {
        double segment_len =
            static_cast<double>(anchors[i + 1].read_pos - anchors[i].read_pos);
        total_bases += segment_len;

        bool cur_in_roi = roi.count(anchors[i].node_id) > 0;
        bool next_in_roi = roi.count(anchors[i + 1].node_id) > 0;
        if (cur_in_roi && next_in_roi) {
          roi_bases += segment_len;
        }
      }

      result.roi_overlap = (total_bases > 0.0) ? roi_bases / total_bases : 0.0;
      above_threshold = result.roi_overlap >= config_.roi_overlap_threshold;
    }

    if (config_.classify_mode == "enrich") {
      result.roi_keep = above_threshold;
    } else {  // deplete
      result.roi_keep = !above_threshold;
    }
  }
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

    /* Write to result file if configured */
    if (config_.result_writer) {
      config_.result_writer->write(map_result, read.read_id, read.len_raw_signal);
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

}  // namespace piru::mapping
