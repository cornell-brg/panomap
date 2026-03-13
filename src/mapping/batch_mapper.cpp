// SPDX-License-Identifier: MIT

#include "mapping/batch_mapper.hpp"

#include <iostream>
#include <string>
#include <utility>

#include "mapping/graph_chainer.hpp"
#include "mapping/path_chainer.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

void SeedLookup::lookup(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& out_hits) const {
  if (!store_) return;
  out_hits.clear();
  out_hits.reserve(seeds.seeds.size());
  for (const auto& seed : seeds.seeds) {
    const auto* hits = store_->lookup(seed.hash);
    if (!hits) continue;
    if (hits->size() > freq_threshold_) continue;  // skip overly frequent seeds
    for (const auto& h : *hits) {
      out_hits.push_back(NodeAnchor{
          .node_id = static_cast<std::uint32_t>(h.node_id),
          .offset = static_cast<std::uint32_t>(h.offset),
          .read_pos = static_cast<std::uint32_t>(seed.position),
          .span = static_cast<std::uint16_t>(std::min(seed.length, std::size_t{0xFFFF})),
          .length = static_cast<std::uint16_t>(std::min(h.length, std::size_t{0xFFFF})),
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

  // Create unified event pipeline (event detection + normalization)
  comps.event_pipeline = signal::make_event_pipeline(config_.event_pipeline_config);

  comps.fuzzy_quantizer = signal::make_fuzzy_quantizer(config_.fuzzy_config);
  comps.seed_extractor = signal::make_seed_extractor(config_.seed_config);
  comps.seed_store = config_.seed_store;
  comps.graph_store = config_.graph_store;
  if (!comps.seed_store) {
    throw std::runtime_error("BatchMapper requires a SeedStore for lookup");
  }

  // Create chainer backend
  if (config_.chainer_backend == "path-chain") {
    if (!config_.linearization_coords) {
      throw std::runtime_error("PathChainer requires linearization_coords");
    }
    if (!config_.path_lengths) {
      throw std::runtime_error("PathChainer requires path_lengths for bounds checking");
    }
    auto path_config = PathChainerConfig::from_parsed(config_.chainer_parsed);
    path_config.merge_anchors = config_.enable_anchor_merge;
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
  } else {
    throw std::runtime_error("Unknown chainer backend: " + config_.chainer_backend);
  }
  const std::size_t freq_threshold = comps.seed_store->frequency_threshold();
  comps.lookup = SeedLookup(comps.seed_store, freq_threshold, config_.max_total_hits);

  // Log pipeline configuration
  LOG_DEBUG("Pipeline: " + comps.chainer->name() + " chaining");

  // Create result formatter if we have a graph store (needed for result output)
  const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);
  if (adj_store && config_.result_writer) {
    comps.result_formatter =
        std::make_unique<ResultFormatter>(adj_store->graph(), config_.formatter_config);
    LOG_INFO("Result formatter enabled for output (min_secondary_ratio=" +
             std::to_string(config_.formatter_config.min_secondary_ratio) + ")");
  }

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

    // Accumulate output stats
    auto batch_stats = output_batch(batch);
    stats.reads_mapped += batch_stats.reads_mapped;
    stats.reads_unmapped += batch_stats.reads_unmapped;
    stats.results_written += batch_stats.results_written;
    stats.primary_alignments += batch_stats.primary_alignments;
    stats.secondary_alignments += batch_stats.secondary_alignments;

    batch.clear();
  }

  // Log detailed summary
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
  // Calculate estimated memory before processing
  std::size_t est_mem_mb = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    est_mem_mb += batch.raw_reads[i].len_raw_signal * sizeof(int16_t) / 1024 / 1024;
  }
  LOG_INFO("Processing batch: " + std::to_string(batch.num_reads) + " reads, ~" +
           std::to_string(est_mem_mb) + " MB raw signal");

  executor_->parallel_for(0, batch.num_reads, 1, [&](std::size_t i) { process_read(batch, i); });

  // Calculate seed hits memory
  std::size_t total_hits = 0;
  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    total_hits += batch.seed_hits[i].size();
  }
  const std::size_t hits_mem_mb = total_hits * sizeof(NodeAnchor) / 1024 / 1024;
  LOG_INFO("Batch complete: total_hits=" + std::to_string(total_hits) + " (~" +
           std::to_string(hits_mem_mb) + " MB)");
}

void BatchMapper::process_read(BatchBuffer& batch, std::size_t index) const {
  const auto& read = batch.raw_reads[index];

  // Signal processing: event detection + normalization
  batch.normalized[index] = components_.event_pipeline->process(batch.raw_reads[index]);

  // Debug: log raw signal size vs event count for each read
  // Expected: events ≈ basepairs (each event ~1bp), raw_signal ≈ events * samples_per_base (~9)
  LOG_DEBUG("read=" + read.read_id + " raw_samples=" + std::to_string(read.len_raw_signal) +
            " events=" + std::to_string(batch.normalized[index].samples.size()) + " ratio=" +
            std::to_string(batch.normalized[index].samples.empty()
                               ? 0.0
                               : static_cast<double>(read.len_raw_signal) /
                                     batch.normalized[index].samples.size()));

  batch.fuzzy_quantized[index] = components_.fuzzy_quantizer->quantize(batch.normalized[index]);
  batch.seeds[index] = components_.seed_extractor->extract(batch.fuzzy_quantized[index]);

  // Lookup seeds in the index and collect hits.
  components_.lookup.lookup(batch.seeds[index], batch.seed_hits[index]);

  // ROI anchor filtering: keep only hits whose node_id is in the ROI set
  if (config_.roi_filter_anchors && config_.roi_nodes) {
    const auto& roi = *config_.roi_nodes;
    auto& hits = batch.seed_hits[index];
    hits.erase(std::remove_if(hits.begin(), hits.end(),
                               [&roi](const NodeAnchor& a) {
                                 return roi.count(a.node_id) == 0;
                               }),
               hits.end());
  }

  // Chain: expands NodeAnchors to PathAnchors internally, then DP chains
  ChainResult chain_result = components_.chainer->chain(batch.seed_hits[index]);

  // Free seed hits immediately -- they're not needed after chaining and can be
  // very large (especially with kmer seeds). Without this, all reads in the
  // batch hold their hits simultaneously, causing OOM on large batches.
  batch.seed_hits[index].clear();
  batch.seed_hits[index].shrink_to_fit();

  // Build unified map result from chain result
  ReadMapResult& result = batch.map_results[index];
  result.mappings.clear();
  result.expanded_anchor_count = chain_result.expanded_anchor_count;

  for (const auto& chain : chain_result.chains) {
    Mapping mapping;
    mapping.anchors = chain.anchors;
    mapping.chain_score = chain.score;
    result.mappings.push_back(std::move(mapping));
  }

  // ROI classification (if --roi is active)
  if (config_.roi_nodes && !config_.classify_mode.empty() && result.mapped()) {
    bool above_threshold;

    if (config_.roi_filter_anchors) {
      // --chain-target: classify by chain score threshold.
      // All anchors are ROI nodes, so overlap is meaningless.
      double chain_score = result.mappings[0].chain_score;
      result.roi_overlap = chain_score;  // repurpose field for output
      above_threshold = chain_score >= config_.roi_score_threshold;
    } else {
      // --chain-genome: classify by ROI node overlap fraction
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

  // Try to get path names from graph store
  const index::AdjListGraphStore* adj_store =
      dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);

  for (std::size_t i = 0; i < batch.num_reads; ++i) {
    const auto& read = batch.raw_reads[i];
    const auto& map_result = batch.map_results[i];

    // Track mapped/unmapped
    const bool is_mapped = map_result.mapped();
    if (is_mapped) {
      ++stats.reads_mapped;
    } else {
      ++stats.reads_unmapped;
    }

    // Write to result file if configured
    if (config_.result_writer && components_.result_formatter) {
      auto results =
          components_.result_formatter->format(map_result, read.read_id, read.len_raw_signal);
      for (std::size_t r = 0; r < results.size(); ++r) {
        config_.result_writer->write(results[r]);
        ++stats.results_written;
        if (r == 0) {
          ++stats.primary_alignments;
        } else {
          ++stats.secondary_alignments;
        }
      }
    }

    // Suppress unused variable warnings
    (void)adj_store;
  }

  return stats;
}

}  // namespace piru::mapping
