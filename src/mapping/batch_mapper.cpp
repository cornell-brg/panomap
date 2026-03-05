// SPDX-License-Identifier: MIT

#include "mapping/batch_mapper.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "mapping/anchor_expander.hpp"
#include "mapping/anchor_merger.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Dump anchors to file for visualization (first read only).
void dumpAnchorsToFile(const char* filename, const std::vector<Anchor>& anchors,
                       const std::string& read_id, const index::GraphStore* graph_store) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        LOG_WARN("Failed to open anchor dump file: " + std::string(filename));
        return;
    }

    // Try to get graph for path metadata
    const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(graph_store);
    if (!adj_store) {
        LOG_WARN("Cannot dump anchors: graph_store is not AdjListGraphStore");
        return;
    }

    const auto& graph = adj_store->graph();
    const auto& paths = graph.paths();

    // Write path metadata header
    for (std::size_t path_id = 0; path_id < paths.size(); ++path_id) {
        const auto& path = paths[path_id];
        out << "#PATH\t" << path_id << "\t" << path.name << "\n";
    }

    // Write anchors header
    out << "#ANCHORS\tread_id\tpath_id\tpath_name\tquery_pos\tref_coord\tlength\tnode_id\n";

    // Write anchor data
    for (const auto& anchor : anchors) {
        std::string path_name;
        if (anchor.path_id < paths.size()) {
            path_name = paths[anchor.path_id].name;
        } else {
            path_name = "path" + std::to_string(anchor.path_id);
        }

        out << read_id << "\t" << anchor.path_id << "\t" << path_name << "\t" << anchor.query_pos
            << "\t" << anchor.ref_coord << "\t" << anchor.length << "\t" << anchor.node_id << "\n";
    }

    out.close();
    LOG_INFO("Dumped " + std::to_string(anchors.size()) + " anchors to " + std::string(filename));
}

// Dump chain (DP output) to file for visualization (first read only).
void dumpChainToFile(const char* filename, const ReadMapResult& result, const std::string& read_id,
                     const index::GraphStore* graph_store) {
    if (result.mappings.empty()) {
        LOG_WARN("Cannot dump chain: no mappings for read " + read_id);
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        LOG_WARN("Failed to open chain dump file: " + std::string(filename));
        return;
    }

    // Try to get graph for path metadata
    const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(graph_store);
    if (!adj_store) {
        LOG_WARN("Cannot dump chain: graph_store is not AdjListGraphStore");
        return;
    }

    const auto& graph = adj_store->graph();
    const auto& paths = graph.paths();

    // Use primary mapping
    const auto& primary = result.mappings[0];

    // Write chain metadata header
    out << "#CHAIN\tread_id\tscore\tchained_anchors\n";
    out << "#CHAIN\t" << read_id << "\t" << primary.chain_score << "\t" << primary.anchors.size()
        << "\n";

    // Write path metadata header
    for (std::size_t path_id = 0; path_id < paths.size(); ++path_id) {
        const auto& path = paths[path_id];
        out << "#PATH\t" << path_id << "\t" << path.name << "\n";
    }

    // Write chained anchors header
    out << "#ANCHORS\tread_id\tpath_id\tpath_name\tquery_pos\tref_coord\tlength\tnode_id\tscore\n";

    // Write chained anchor data
    for (const auto& anchor : primary.anchors) {
        std::string path_name;
        if (anchor.path_id < paths.size()) {
            path_name = paths[anchor.path_id].name;
        } else {
            path_name = "path" + std::to_string(anchor.path_id);
        }

        out << read_id << "\t" << anchor.path_id << "\t" << path_name << "\t" << anchor.read_pos
            << "\t" << anchor.ref_coord << "\t" << anchor.target.length << "\t"
            << anchor.target.node_id << "\t" << anchor.score << "\n";
    }

    out.close();
    LOG_INFO("Dumped chain with " + std::to_string(primary.anchors.size()) + " anchors to " +
             std::string(filename));
}

// Dump ALL seeds for a read (including those with no hits in the index).
// This is critical for debugging: we need to see which hashes the read produces
// and whether they exist in the index at all.
void dumpAllReadSeedsToFile(const char* filename, const std::string& read_id,
                            const signal::SeedBuffer& seeds, const index::SeedStore* seed_store,
                            std::size_t freq_threshold) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        LOG_WARN("Failed to open read seeds file: " + std::string(filename));
        return;
    }

    // Write header
    out << "#READ_SEEDS\tread_id\tnum_seeds\tfreq_threshold\n";
    out << "#READ_SEEDS\t" << read_id << "\t" << seeds.seeds.size() << "\t" << freq_threshold
        << "\n";

    // Write per-seed details
    out << "hash\tread_pos\tlength\tin_index\tfrequency\tfiltered\n";

    for (const auto& seed : seeds.seeds) {
        const auto* hits = seed_store ? seed_store->lookup(seed.hash) : nullptr;
        const bool in_index = (hits != nullptr);
        const std::size_t frequency = in_index ? hits->size() : 0;
        const bool filtered = in_index && (frequency > freq_threshold);

        out << std::hex << seed.hash << std::dec << "\t" << seed.position << "\t" << seed.length
            << "\t" << (in_index ? "yes" : "no") << "\t" << frequency << "\t"
            << (filtered ? "yes" : "no") << "\n";
    }

    out.close();
}

// Dump hit statistics for a read to analyze frequency distributions.
void dumpHitStatsToFile(const char* filename, const std::string& read_id,
                        const signal::SeedBuffer& seeds, const std::vector<SeedHitRecord>& hits,
                        std::size_t freq_threshold) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        LOG_WARN("Failed to open hit stats file: " + std::string(filename));
        return;
    }

    // Collect frequency statistics
    std::vector<std::size_t> frequencies;
    frequencies.reserve(hits.size());
    for (const auto& hit : hits) {
        frequencies.push_back(hit.frequency);
    }

    // Sort for percentiles
    std::sort(frequencies.begin(), frequencies.end());

    // Compute stats
    std::size_t total_hits = hits.size();
    std::size_t min_freq = frequencies.empty() ? 0 : frequencies.front();
    std::size_t max_freq = frequencies.empty() ? 0 : frequencies.back();
    std::size_t median_freq = frequencies.empty() ? 0 : frequencies[frequencies.size() / 2];
    std::size_t p90_freq = frequencies.empty() ? 0 : frequencies[frequencies.size() * 9 / 10];

    // Count unique hashes that got hits vs total seeds
    std::unordered_set<std::uint64_t> hashes_with_hits;
    for (const auto& hit : hits) {
        hashes_with_hits.insert(hit.hash);
    }

    // Write summary header
    out << "#SUMMARY\tread_id\tnum_seeds\tseeds_with_hits\ttotal_hits\tfreq_threshold\n";
    out << "#SUMMARY\t" << read_id << "\t" << seeds.seeds.size() << "\t" << hashes_with_hits.size()
        << "\t" << total_hits << "\t" << freq_threshold << "\n";

    out << "#FREQ_STATS\tmin\tmax\tmedian\tp90\n";
    out << "#FREQ_STATS\t" << min_freq << "\t" << max_freq << "\t" << median_freq << "\t"
        << p90_freq << "\n";

    // Write per-seed details header
    out << "#SEEDS\thash\tread_pos\thit_count\n";

    // Group hits by hash to get per-seed hit counts
    std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>>
        hash_info;  // hash -> (read_pos, freq)
    for (const auto& hit : hits) {
        if (hash_info.find(hit.hash) == hash_info.end()) {
            hash_info[hit.hash] = {hit.read_pos, hit.frequency};
        }
    }

    // Write per-seed info
    for (const auto& [hash, info] : hash_info) {
        out << std::hex << hash << std::dec << "\t" << info.first << "\t" << info.second << "\n";
    }

    out.close();
}

}  // namespace

void SeedLookup::lookup(const signal::SeedBuffer& seeds,
                        std::vector<SeedHitRecord>& out_hits) const {
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
                .chain_id =
                    graph_store_ ? graph_store_->chainId(h.node_id) : std::optional<std::int64_t>{},
                .linear_pos = graph_store_ ? graph_store_->linearPosition(h.node_id)
                                           : std::optional<std::int64_t>{},
                .frequency = hits->size(),
            });
        }
    }
}

void BatchBuffer::resize(std::size_t capacity) {
    raw_reads.resize(capacity);
    normalized.resize(capacity);
    fuzzy_quantized.resize(capacity);
    seeds.resize(capacity);
    seed_hits.resize(capacity);
    anchors.resize(capacity);
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
        anchors[i].clear();
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

    // Create AnchorExpander based on linearization type
    if (config_.linearization_coords) {
        // Path-walk linearization -> PathWalkExpander
        if (!config_.path_lengths) {
            throw std::runtime_error("PathWalkExpander requires path_lengths for bounds checking");
        }
        comps.expander = std::make_unique<PathWalkExpander>(*config_.linearization_coords,
                                                            *config_.path_lengths);
    } else {
        // Superbubble linearization -> SuperbubbleExpander
        if (!config_.graph_store) {
            throw std::runtime_error("BatchMapper requires GraphStore for superbubble expansion");
        }
        comps.expander = std::make_unique<SuperbubbleExpander>(config_.graph_store);
    }

    comps.chainer = make_chainer(config_.chainer_backend, config_.chainer_parsed);
    const std::size_t freq_threshold = comps.seed_store->frequency_threshold();
    // Limit the lookup helper to what the SeedStore exposes.
    comps.lookup = SeedLookup(comps.seed_store, comps.graph_store, freq_threshold);

    // Log pipeline configuration
    LOG_DEBUG("Pipeline: " + comps.expander->name() + " expansion + " + comps.chainer->name() +
              " chaining");

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

    LOG_DEBUG("BatchMapper starting: batch_capacity_reads=" +
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
    const std::size_t hits_mem_mb = total_hits * sizeof(SeedHitRecord) / 1024 / 1024;
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

    // Dump hit stats if --dump-hit-stats is set
    static std::atomic<std::size_t> hit_stats_dump_counter{0};
    if (!config_.dump_hit_stats_dir.empty()) {
        std::size_t dump_num = hit_stats_dump_counter.fetch_add(1);
        std::string hit_stats_file =
            config_.dump_hit_stats_dir + "/read_" + std::to_string(dump_num) + "_hits.tsv";
        std::size_t freq_threshold =
            config_.seed_store ? config_.seed_store->frequency_threshold() : 0;
        dumpHitStatsToFile(hit_stats_file.c_str(), read.read_id, batch.seeds[index],
                           batch.seed_hits[index], freq_threshold);
    }

    // Dump ALL read seeds (including no-hit) if --dump-read-seeds is set
    static std::atomic<std::size_t> read_seeds_dump_counter{0};
    if (!config_.dump_read_seeds_dir.empty()) {
        std::size_t dump_num = read_seeds_dump_counter.fetch_add(1);
        std::string seeds_file =
            config_.dump_read_seeds_dir + "/read_" + std::to_string(dump_num) + "_seeds.tsv";
        std::size_t freq_threshold =
            config_.seed_store ? config_.seed_store->frequency_threshold() : 0;
        dumpAllReadSeedsToFile(seeds_file.c_str(), read.read_id, batch.seeds[index],
                               config_.seed_store, freq_threshold);
    }

    // Expand seed hits to anchors (explicit expansion stage)
    auto anchors = components_.expander->expand(batch.seed_hits[index]);

    // Merge adjacent/overlapping anchors on same path (optional optimization)
    if (config_.enable_anchor_merge) {
        const auto pre_merge_count = anchors.size();
        anchors = AnchorMerger::merge(anchors, AnchorMergerConfig{});
        LOG_DEBUG("Anchor merge: " + std::to_string(pre_merge_count) + " -> " +
                  std::to_string(anchors.size()) + " (" +
                  std::to_string(pre_merge_count > 0
                                     ? 100 * (pre_merge_count - anchors.size()) / pre_merge_count
                                     : 0) +
                  "% reduction)");
    }

    // Dump anchors if --dump-anchors is set
    static std::atomic<std::size_t> anchor_dump_counter{0};
    if (!config_.dump_anchors_dir.empty()) {
        std::size_t dump_num = anchor_dump_counter.fetch_add(1);
        std::string anchor_file =
            config_.dump_anchors_dir + "/read_" + std::to_string(dump_num) + "_anchors.tsv";
        dumpAnchorsToFile(anchor_file.c_str(), anchors, read.read_id, config_.graph_store);
    }

    // Chain anchors
    ChainResult chain_result = components_.chainer->chain(anchors);

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
        const auto& anchors = result.mappings[0].anchors;
        const auto& roi = *config_.roi_nodes;

        double roi_bases = 0.0;
        double total_bases = 0.0;

        for (std::size_t i = 0; i + 1 < anchors.size(); ++i) {
            double segment_len = static_cast<double>(anchors[i + 1].read_pos - anchors[i].read_pos);
            total_bases += segment_len;

            bool cur_in_roi = roi.count(anchors[i].target.node_id) > 0;
            bool next_in_roi = roi.count(anchors[i + 1].target.node_id) > 0;
            if (cur_in_roi && next_in_roi) {
                roi_bases += segment_len;
            }
        }

        result.roi_overlap = (total_bases > 0.0) ? roi_bases / total_bases : 0.0;

        bool above_threshold = result.roi_overlap >= config_.roi_threshold;
        if (config_.classify_mode == "enrich") {
            result.roi_keep = above_threshold;
        } else {  // deplete
            result.roi_keep = !above_threshold;
        }
    }

    // Dump chain if --dump-chains is set
    static std::atomic<std::size_t> chain_dump_counter{0};
    if (!config_.dump_chains_dir.empty()) {
        std::size_t dump_num = chain_dump_counter.fetch_add(1);
        std::string chain_file =
            config_.dump_chains_dir + "/read_" + std::to_string(dump_num) + "_chain.tsv";
        dumpChainToFile(chain_file.c_str(), result, read.read_id, config_.graph_store);
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
