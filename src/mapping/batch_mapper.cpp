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
#include <utility>

#include "alignment/signal_utils.hpp"
#include "index/signal_store.hpp"
#include "mapping/anchor_expander.hpp"
#include "mapping/anchor_merger.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Compute path length by summing signal sizes for each node on the path.
std::int64_t computePathLength(const index::AlnGraph& graph,
                               const index::AlnPath& path,
                               const index::SignalStore* signal_store) {
    // Build label → node index mapping
    std::unordered_map<std::string, std::size_t> label_to_idx;
    for (std::size_t i = 0; i < graph.nodeCount(); ++i) {
        label_to_idx[graph.node(i).label] = i;
    }

    std::int64_t cumulative = 0;
    for (const auto& step : path.steps) {
        auto it = label_to_idx.find(step.node_id);
        if (it == label_to_idx.end()) continue;

        const std::size_t node_idx = it->second;

        // Use signal size directly (no overlap math needed)
        if (signal_store) {
            const auto* sig = signal_store->get(node_idx);
            if (sig) {
                cumulative += static_cast<std::int64_t>(alignment::signalLength(*sig));
            }
        }
    }
    return cumulative;
}

// Dump anchors to file for visualization (first read only).
void dumpAnchorsToFile(const char* filename,
                       const std::vector<Anchor>& anchors,
                       const std::string& read_id,
                       const index::GraphStore* graph_store,
                       const index::SignalStore* signal_store) {
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
        std::int64_t length = computePathLength(graph, path, signal_store);
        out << "#PATH\t" << path_id << "\t" << path.name << "\t" << length << "\n";
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

        out << read_id << "\t"
            << anchor.path_id << "\t"
            << path_name << "\t"
            << anchor.query_pos << "\t"
            << anchor.ref_coord << "\t"
            << anchor.length << "\t"
            << anchor.node_id << "\n";
    }

    out.close();
    LOG_INFO("Dumped " + std::to_string(anchors.size()) + " anchors to " + std::string(filename));
}

// Dump chain (DP output) to file for visualization (first read only).
void dumpChainToFile(const char* filename,
                     const ReadMapResult& result,
                     const std::string& read_id,
                     const index::GraphStore* graph_store,
                     const index::SignalStore* signal_store) {
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
    out << "#CHAIN\t" << read_id << "\t" << primary.chain_score << "\t" << primary.anchors.size() << "\n";

    // Write path metadata header
    for (std::size_t path_id = 0; path_id < paths.size(); ++path_id) {
        const auto& path = paths[path_id];
        std::int64_t length = computePathLength(graph, path, signal_store);
        out << "#PATH\t" << path_id << "\t" << path.name << "\t" << length << "\n";
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

        out << read_id << "\t"
            << anchor.path_id << "\t"
            << path_name << "\t"
            << anchor.read_pos << "\t"
            << anchor.ref_coord << "\t"
            << anchor.target.length << "\t"
            << anchor.target.node_id << "\t"
            << anchor.score << "\n";
    }

    out.close();
    LOG_INFO("Dumped chain with " + std::to_string(primary.anchors.size()) + " anchors to " + std::string(filename));
}

}  // namespace

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
    normalized.resize(capacity);
    fuzzy_quantized.resize(capacity);
    alignment_quantized.resize(capacity);
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
        alignment_quantized[i] = signal::AlignmentQuantizedSignal{};
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
    LOG_DEBUG("Pipeline: " + comps.expander->name() + " expansion + " +
             clusterer_name + " clustering");

    // Create result formatter if we have a graph store (needed for result output)
    const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);
    if (adj_store && config_.result_writer) {
        comps.result_formatter = std::make_unique<ResultFormatter>(
            adj_store->graph(), config_.signal_store, config_.formatter_config);
        LOG_INFO("Result formatter enabled for output (min_secondary_ratio=" +
                 std::to_string(config_.formatter_config.min_secondary_ratio) + ")");
    }

    // Create chain aligner if alignment is enabled
    if (config_.enable_alignment) {
        comps.chain_aligner = std::make_unique<alignment::ChainAligner>(config_.align_config);
        comps.signal_store = config_.signal_store;
        std::string backend_name =
            config_.align_config.backend == alignment::AlignerBackend::kPathGuided ? "path-guided" :
            config_.align_config.backend == alignment::AlignerBackend::kRadius ? "radius" : "auto";
        LOG_INFO("Signal-level alignment enabled: backend=" + backend_name);
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
    const auto& read = batch.raw_reads[index];

#ifdef PIRU_DUMP_GRAPHS
    // Dump read signal processing stages to file (first read only)
    static std::atomic<int> dump_count{0};
    const int dump_idx = dump_count.fetch_add(1);
    const bool should_dump = (dump_idx < 1);  // Only dump first read

    std::ofstream signal_file;
    if (should_dump) {
        std::string filename = "signal_dump_" + std::to_string(dump_idx) + ".txt";
        signal_file.open(filename);
        if (!signal_file.is_open()) {
            LOG_WARN("Failed to open signal dump file: " + filename);
        }
    }
#endif

    // Signal processing: event detection + normalization
    batch.normalized[index] = components_.event_pipeline->process(batch.raw_reads[index]);

    // Debug: log raw signal size vs event count for each read
    // Expected: events ≈ basepairs (each event ~1bp), raw_signal ≈ events * samples_per_base (~9)
    LOG_DEBUG("read=" + read.read_id +
              " raw_samples=" + std::to_string(read.len_raw_signal) +
              " events=" + std::to_string(batch.normalized[index].samples.size()) +
              " ratio=" + std::to_string(
                  batch.normalized[index].samples.empty() ? 0.0 :
                  static_cast<double>(read.len_raw_signal) / batch.normalized[index].samples.size()));

    batch.fuzzy_quantized[index] = components_.fuzzy_quantizer->quantize(batch.normalized[index]);
    batch.alignment_quantized[index] = components_.alignment_quantizer->quantize(batch.normalized[index]);
    batch.seeds[index] = components_.seed_extractor->extract(batch.fuzzy_quantized[index]);

#ifdef PIRU_DUMP_GRAPHS
    if (should_dump && signal_file.is_open()) {
        // Line 1: Metadata (read name)
        signal_file << read.read_id << "\n";

        // Line 2: Raw signal (ADC values converted to picoamps)
        const auto& raw = read.raw_signal;
        const float digitisation = read.digitisation;
        const float offset = read.offset;
        const float range = read.range;

        for (std::size_t i = 0; i < read.len_raw_signal; ++i) {
            if (i > 0) signal_file << ",";
            float pA = ((raw[i] + offset) * range) / digitisation;
            signal_file << pA;
        }
        signal_file << "\n";

        // Line 3: Normalized signal (output of event pipeline)
        const auto& normalized = batch.normalized[index].samples;
        for (std::size_t i = 0; i < normalized.size(); ++i) {
            if (i > 0) signal_file << ",";
            signal_file << normalized[i];
        }
        signal_file << "\n";

        // Line 4: Fuzzy quantized tokens
        const auto& fuzzy = batch.fuzzy_quantized[index].tokens;
        for (std::size_t i = 0; i < fuzzy.size(); ++i) {
            if (i > 0) signal_file << ",";
            signal_file << static_cast<int>(fuzzy[i]);
        }
        signal_file << "\n";

        signal_file.close();
    }
#endif

    // Lookup seeds in the index and collect hits.
    components_.lookup.lookup(batch.seeds[index], batch.seed_hits[index]);

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
        std::string anchor_file = config_.dump_anchors_dir + "/read_" + std::to_string(dump_num) + "_anchors.tsv";
        dumpAnchorsToFile(anchor_file.c_str(), anchors, read.read_id, config_.graph_store, config_.signal_store);
    }

    // Cluster anchors
    ClusterSummary summary = components_.clusterer->cluster(anchors);

    // Build unified map result from cluster summary
    ReadMapResult& result = batch.map_results[index];
    result.mappings.clear();
    result.expanded_anchor_count = summary.expanded_anchor_count;

    for (const auto& cluster : summary.clusters) {
        Mapping mapping;
        mapping.anchors = cluster.anchors;
        mapping.chain_score = cluster.cluster_score;

        // Run alignment if enabled
        if (components_.chain_aligner && components_.signal_store && config_.graph_store &&
            cluster.anchors.size() >= 2) {
            // Convert SeedAnchors to alignment::Anchors
            std::vector<alignment::Anchor> align_anchors;
            align_anchors.reserve(cluster.anchors.size());
            for (const auto& seed_anchor : cluster.anchors) {
                alignment::Anchor a;
                a.graph_pos.node_id = static_cast<std::uint32_t>(seed_anchor.target.node_id);
                a.graph_pos.offset = static_cast<std::uint32_t>(seed_anchor.target.offset);
                a.query_pos = static_cast<std::uint32_t>(seed_anchor.read_pos);
                align_anchors.push_back(a);
            }

            // Run alignment
            auto align_result = components_.chain_aligner->align(
                *config_.graph_store, *components_.signal_store,
                batch.alignment_quantized[index], align_anchors);

            if (align_result.valid()) {
                mapping.alignment_cost = align_result.total_cost;
                std::size_t query_len = cluster.anchors.back().read_pos - cluster.anchors.front().read_pos;
                if (query_len > 0) {
                    mapping.normalized_alignment_cost = align_result.normalizedCost(query_len);
                }
                mapping.alignment_path = std::move(align_result.path);
                mapping.segments_aligned = align_result.segments_aligned;
            }
        }

        result.mappings.push_back(std::move(mapping));
    }

    // Dump chain if --dump-chains is set
    static std::atomic<std::size_t> chain_dump_counter{0};
    if (!config_.dump_chains_dir.empty()) {
        std::size_t dump_num = chain_dump_counter.fetch_add(1);
        std::string chain_file = config_.dump_chains_dir + "/read_" + std::to_string(dump_num) + "_chain.tsv";
        dumpChainToFile(chain_file.c_str(), result, read.read_id, config_.graph_store, config_.signal_store);
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
            auto results = components_.result_formatter->format(
                map_result, read.read_id, read.len_raw_signal);
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

        // Debug output to stdout (disabled - use -o for PAF/GAF output)
        // output_ << read.read_id
        //         << "\tseeds=" << seeds_for_read.size()
        //         << "\thits=" << hits_for_read.size()
        //         << "\tanchors=" << clusters_for_read.expanded_anchor_count
        //         << "\tchains=" << clusters_for_read.clusters.size()
        //         << "\tlen=" << read.len_raw_signal
        //         << "\n";
        //
        // // Output all chains (one per line)
        // for (std::size_t chain_idx = 0; chain_idx < clusters_for_read.clusters.size(); ++chain_idx) {
        //     const auto& chain = clusters_for_read.clusters[chain_idx];
        //     if (chain.anchors.empty()) continue;
        //
        //     const auto& first_anchor = chain.anchors.front();
        //     const auto& last_anchor = chain.anchors.back();
        //
        //     // Get path name
        //     std::string path_name;
        //     if (adj_store && first_anchor.path_id < adj_store->graph().pathCount()) {
        //         path_name = adj_store->graph().paths()[first_anchor.path_id].name;
        //     } else {
        //         path_name = "path" + std::to_string(first_anchor.path_id);
        //     }
        //
        //     // Compute intervals (start of first anchor to end of last anchor)
        //     const auto q_start = first_anchor.read_pos;
        //     const auto q_end = last_anchor.read_pos + last_anchor.target.length;
        //     const auto r_start = first_anchor.ref_coord;
        //     const auto r_end = last_anchor.ref_coord + last_anchor.target.length;
        //
        //     output_ << "  chain" << chain_idx << "\t"
        //             << path_name << "\t"
        //             << "score=" << static_cast<int>(chain.cluster_score) << "\t"
        //             << "n=" << chain.anchors.size() << "\t"
        //             << "q[" << q_start << "-" << q_end << "]\t"
        //             << "r[" << r_start << "-" << r_end << "]\n";
        // }
        //
        // // Blank line between reads
        // output_ << "\n";

        // Suppress unused variable warnings
        (void)adj_store;
    }

    return stats;
}

}  // namespace piru::mapping
