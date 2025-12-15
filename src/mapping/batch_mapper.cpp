// SPDX-License-Identifier: MIT

#include "mapping/batch_mapper.hpp"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>

#include "mapping/anchor_expander.hpp"
#include "mapping/anchor_merger.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

void SeedLookup::lookup(const signal::SeedBuffer& seeds, std::vector<SeedHitRecord>& out_hits) const {
    if (!store_) return;
    out_hits.clear();
    out_hits.reserve(seeds.seeds.size());
    for (const auto& seed : seeds.seeds) {
        const auto* hits = store_->lookup(seed.hash);
        if (!hits) continue;
        if (hits->size() > freq_threshold_) continue;  // skip overly frequent seeds
        for (const auto& h : *hits) {
            out_hits.push_back(SeedHitRecord{
                .target = h,
                .read_pos = seed.position,
                .hash = seed.hash,
                .span = seed.length,  // Use length (may differ from span after merging)
                .chain_id = graph_store_ ? graph_store_->chainId(h.node_id) : std::optional<std::int64_t>{},
                .linear_pos = graph_store_ ? graph_store_->linearPosition(h.node_id) : std::optional<std::int64_t>{},
                .frequency = hits->size(),
            });
        }
    }
}

void BatchBuffer::resize(std::size_t capacity) {
    raw_reads.resize(capacity);
    events.resize(capacity);
    normalized.resize(capacity);
    fuzzy_quantized.resize(capacity);
    alignment_quantized.resize(capacity);
    seeds.resize(capacity);
    seed_hits.resize(capacity);
    anchors.resize(capacity);
    clusters.resize(capacity);
    alignment_notes.resize(capacity);
    num_reads = 0;
}

void BatchBuffer::clear() {
    for (std::size_t i = 0; i < num_reads; ++i) {
        raw_reads[i] = io::RawRead{};
        events[i].events.clear();
        normalized[i] = signal::NormalizedSignal{};
        fuzzy_quantized[i].tokens.clear();
        alignment_quantized[i] = signal::AlignmentQuantizedSignal{};
        seeds[i].seeds.clear();
        seed_hits[i].clear();
        anchors[i].clear();
        clusters[i] = ClusterSummary{};
        alignment_notes[i].clear();
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
    comps.event_detector = signal::make_event_detector(config_.event_config);
    comps.normalizer = signal::make_signal_normalizer(config_.normalizer_config);
    comps.fuzzy_quantizer = signal::make_fuzzy_quantizer(config_.fuzzy_config);
    comps.alignment_quantizer = signal::make_alignment_quantizer(config_.alignment_config);
    comps.seed_extractor = signal::make_seed_extractor(config_.seed_config);
    comps.seed_store = config_.seed_store;
    comps.graph_store = config_.graph_store;
    if (!comps.seed_store) {
        throw std::runtime_error("BatchMapper requires a SeedStore for lookup");
    }

    // Create AnchorExpander based on linearization type
    if (config_.linearization_coords) {
        // Path-walk linearization → PathWalkExpander
        comps.expander = std::make_unique<PathWalkExpander>(*config_.linearization_coords);
    } else {
        // Superbubble linearization → SuperbubbleExpander
        if (!config_.graph_store) {
            throw std::runtime_error("BatchMapper requires GraphStore for superbubble expansion");
        }
        comps.expander = std::make_unique<SuperbubbleExpander>(config_.graph_store);
    }

    comps.clusterer = make_anchor_clusterer(
        config_.clusterer_config,
        config_.graph_store);
    const std::size_t freq_threshold = comps.seed_store->frequency_threshold();
    // Limit the lookup helper to what the SeedStore exposes.
    comps.lookup = SeedLookup(comps.seed_store, comps.graph_store, freq_threshold);

    // Pipeline validation warnings
    const bool is_path_walk = (config_.linearization_coords != nullptr);
    const std::string clusterer_name = comps.clusterer->name();
    const bool is_dp_chain = (clusterer_name == "dp-chain");

    if (is_path_walk && !is_dp_chain) {
        LOG_WARN("Using path-walk linearization with '" + clusterer_name +
                 "' clusterer. Consider --clusterer dp-chain for optimal path-walk support.");
    }
    if (!is_path_walk && is_dp_chain) {
        LOG_WARN("Using DP chaining with superbubble linearization. " +
                 std::string("DP chaining works best with path-walk linearization."));
    }

    // Log pipeline configuration
    LOG_INFO("Pipeline: " + comps.expander->name() + " expansion + " +
             clusterer_name + " clustering");

    return comps;
}

void BatchMapper::run_alignment_stub(const ClusterSummary& summary, const io::RawRead& read,
                                     std::string& note) const {
    if (summary.anchors.empty()) {
        note = "align=none";
        return;
    }
    const auto backend = components_.clusterer ? components_.clusterer->name() : "unknown";
    note = "align=" + backend + " chains=" + std::to_string(summary.anchors.size());
    if (!read.read_id.empty()) {
        note += " read=" + read.read_id;
    }
}

BatchMapperStats BatchMapper::process_all() {
    BatchMapperStats stats;
    BatchBuffer batch;
    batch.resize(config_.batch_capacity_reads);

    LOG_INFO("BatchMapper starting: batch_capacity_reads=" +
             std::to_string(config_.batch_capacity_reads) +
             ", threads=" + std::to_string(executor_->max_concurrency()) +
             ", fuzzy=" + components_.fuzzy_quantizer->name() +
             ", seeds=" + components_.seed_extractor->name());

    while (true) {
        load_batch(batch);
        if (batch.num_reads == 0) break;

        ++stats.batches;
        stats.reads_processed += batch.num_reads;
        process_batch(batch);
        output_batch(batch);
        batch.clear();
    }

    LOG_INFO("BatchMapper finished: batches=" + std::to_string(stats.batches) +
             ", reads=" + std::to_string(stats.reads_processed));
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
    LOG_INFO("Processing batch: " + std::to_string(batch.num_reads) + " reads, ~" + std::to_string(est_mem_mb) + " MB raw signal");

    executor_->parallel_for(0, batch.num_reads, 1, [&](std::size_t i) { process_read(batch, i); });

    // Calculate seed hits memory
    std::size_t total_hits = 0;
    for (std::size_t i = 0; i < batch.num_reads; ++i) {
        total_hits += batch.seed_hits[i].size();
    }
    const std::size_t hits_mem_mb = total_hits * sizeof(SeedHitRecord) / 1024 / 1024;
    LOG_INFO("Batch complete: total_hits=" + std::to_string(total_hits) + " (~" + std::to_string(hits_mem_mb) + " MB)");
}

void BatchMapper::process_read(BatchBuffer& batch, std::size_t index) const {
    batch.events[index] = components_.event_detector->detect(batch.raw_reads[index]);
    batch.normalized[index] = components_.normalizer->normalize(batch.raw_reads[index], &batch.events[index]);

    batch.fuzzy_quantized[index] = components_.fuzzy_quantizer->quantize(batch.normalized[index], &batch.events[index]);
    batch.alignment_quantized[index] = components_.alignment_quantizer->quantize(batch.normalized[index], &batch.events[index]);
    batch.seeds[index] = components_.seed_extractor->extract(batch.fuzzy_quantized[index], &batch.events[index]);

    // Lookup seeds in the index and collect hits.
    components_.lookup.lookup(batch.seeds[index], batch.seed_hits[index]);

    // Expand seed hits to anchors (explicit expansion stage)
    auto anchors = components_.expander->expand(batch.seed_hits[index]);

    // Merge adjacent/overlapping anchors on same path (optional optimization)
    anchors = AnchorMerger::merge(anchors, AnchorMergerConfig{});

    // Store anchors for debugging (optional field)
    batch.anchors[index] = anchors;

    // Cluster anchors
    batch.clusters[index] = components_.clusterer->cluster(anchors);

    // Alignment stub: record what would be aligned (backend + anchor count).
    run_alignment_stub(batch.clusters[index], batch.raw_reads[index], batch.alignment_notes[index]);
}

void BatchMapper::output_batch(const BatchBuffer& batch) const {
    for (std::size_t i = 0; i < batch.num_reads; ++i) {
        const auto& read = batch.raw_reads[i];
        const auto& seeds_for_read = batch.seeds[i].seeds;
        const auto& hits_for_read = batch.seed_hits[i];
        const auto& clusters_for_read = batch.clusters[i];
        output_ << read.read_id
                << "\tseeds=" << seeds_for_read.size()
                << "\thits=" << hits_for_read.size()
                << "\tanchors=" << clusters_for_read.expanded_anchor_count
                << "\tchains=" << clusters_for_read.anchors.size()
                << "\tlen=" << read.len_raw_signal;

        // Output chain path (node IDs)
        if (!clusters_for_read.anchors.empty()) {
            output_ << "\tpath=";
            for (std::size_t j = 0; j < clusters_for_read.anchors.size(); ++j) {
                if (j > 0) output_ << ",";
                output_ << clusters_for_read.anchors[j].target.node_id;
            }
        }

        output_ << "\n";
    }
}

}  // namespace piru::mapping
