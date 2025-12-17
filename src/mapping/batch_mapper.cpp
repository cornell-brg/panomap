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

#include "mapping/anchor_expander.hpp"
#include "mapping/anchor_merger.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Compute path length by walking path steps and summing (node_length - overlap).
std::int64_t computePathLength(const index::AlnGraph& graph, const index::AlnPath& path) {
    // Build label → node index mapping
    std::unordered_map<std::string, std::size_t> label_to_idx;
    for (std::size_t i = 0; i < graph.nodeCount(); ++i) {
        label_to_idx[graph.node(i).label] = i;
    }

    std::int64_t cumulative = 0;
    for (std::size_t step_idx = 0; step_idx < path.steps.size(); ++step_idx) {
        const auto& step = path.steps[step_idx];
        auto it = label_to_idx.find(step.node_id);
        if (it == label_to_idx.end()) continue;

        const auto& node = graph.node(it->second);
        const std::size_t node_len = node.sequence.size();
        std::size_t overlap = 0;

        if (step_idx < path.steps.size() - 1 && step_idx < path.overlaps.size()) {
            overlap = path.overlaps[step_idx];
        }

        cumulative += static_cast<std::int64_t>(node_len) - static_cast<std::int64_t>(overlap);
    }
    return cumulative;
}

// Dump anchors to file for visualization (first read only).
void dumpAnchorsToFile(const char* filename,
                       const std::vector<Anchor>& anchors,
                       const std::string& read_id,
                       const index::GraphStore* graph_store) {
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
        std::int64_t length = computePathLength(graph, path);
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
                     const ClusterSummary& chain,
                     const std::string& read_id,
                     const index::GraphStore* graph_store) {
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

    // Write chain metadata header
    out << "#CHAIN\tread_id\tscore\tchained_anchors\n";
    out << "#CHAIN\t" << read_id << "\t" << chain.score << "\t" << chain.anchors.size() << "\n";

    // Write path metadata header
    for (std::size_t path_id = 0; path_id < paths.size(); ++path_id) {
        const auto& path = paths[path_id];
        std::int64_t length = computePathLength(graph, path);
        out << "#PATH\t" << path_id << "\t" << path.name << "\t" << length << "\n";
    }

    // Write chained anchors header
    out << "#ANCHORS\tread_id\tpath_id\tpath_name\tquery_pos\tref_coord\tlength\tnode_id\tscore\n";

    // Write chained anchor data
    for (const auto& anchor : chain.anchors) {
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
    LOG_INFO("Dumped chain with " + std::to_string(chain.anchors.size()) + " anchors to " + std::string(filename));
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

    // Create result converter if we have a graph store (needed for result output)
    const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);
    if (adj_store && config_.result_writer) {
        comps.result_converter = std::make_unique<ChainResultConverter>(adj_store->graph());
        LOG_INFO("Result converter enabled for output");
    }

    return comps;
}

void BatchMapper::run_alignment_stub(const ClusterSummary& summary, const io::RawRead& read,
                                     std::string& note) const {
    if (summary.anchors.empty()) {
        note = "align=none";
        return;
    }
    const auto backend = components_.clusterer ? components_.clusterer->name() : "unknown";
    note = "align=" + backend + " chained_anchors=" + std::to_string(summary.anchors.size());
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

    batch.events[index] = components_.event_detector->detect(batch.raw_reads[index]);
    batch.normalized[index] = components_.normalizer->normalize(batch.events[index]);
    batch.fuzzy_quantized[index] = components_.fuzzy_quantizer->quantize(batch.normalized[index]);
    batch.alignment_quantized[index] = components_.alignment_quantizer->quantize(batch.normalized[index]);
    batch.seeds[index] = components_.seed_extractor->extract(batch.fuzzy_quantized[index]);

#ifdef PIRU_DUMP_GRAPHS
    if (should_dump && signal_file.is_open()) {
        // Line 1: Metadata (read name)
        signal_file << read.read_id << "\n";

        // Line 2: Raw signal (ADC values converted to picoamps)
        // Convert ADC to picoamps the same way event detector does
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

        // Line 3: Event signal (event means in picoamps)
        const auto& events = batch.events[index].events;
        for (std::size_t i = 0; i < events.size(); ++i) {
            if (i > 0) signal_file << ",";
            signal_file << events[i].mean;
        }
        signal_file << "\n";

        // Line 4: Normalized signal
        const auto& normalized = batch.normalized[index].samples;
        for (std::size_t i = 0; i < normalized.size(); ++i) {
            if (i > 0) signal_file << ",";
            signal_file << normalized[i];
        }
        signal_file << "\n";

        // Line 5: Fuzzy quantized tokens
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
    anchors = AnchorMerger::merge(anchors, AnchorMergerConfig{});

    // Dump anchors for first read if PIRU_DUMP_ANCHORS is set
    static std::atomic<bool> anchor_dump_done{false};
    const char* anchor_dump_file = std::getenv("PIRU_DUMP_ANCHORS");
    if (anchor_dump_file && !anchor_dump_done.exchange(true)) {
        dumpAnchorsToFile(anchor_dump_file, anchors, read.read_id, config_.graph_store);
    }

    // Cluster anchors
    batch.clusters[index] = components_.clusterer->cluster(anchors);

    // Dump chain for first read if PIRU_DUMP_CHAIN is set
    static std::atomic<bool> chain_dump_done{false};
    const char* chain_dump_file = std::getenv("PIRU_DUMP_CHAIN");
    if (chain_dump_file && !chain_dump_done.exchange(true)) {
        dumpChainToFile(chain_dump_file, batch.clusters[index], read.read_id, config_.graph_store);
    }

    // Alignment stub: record what would be aligned (backend + anchor count).
    run_alignment_stub(batch.clusters[index], batch.raw_reads[index], batch.alignment_notes[index]);
}

void BatchMapper::output_batch(const BatchBuffer& batch) const {
    // Try to get path names from graph store
    const index::AdjListGraphStore* adj_store =
        dynamic_cast<const index::AdjListGraphStore*>(config_.graph_store);

    for (std::size_t i = 0; i < batch.num_reads; ++i) {
        const auto& read = batch.raw_reads[i];
        const auto& seeds_for_read = batch.seeds[i].seeds;
        const auto& hits_for_read = batch.seed_hits[i];
        const auto& clusters_for_read = batch.clusters[i];

        // Write to result file if configured
        if (config_.result_writer && components_.result_converter) {
            auto results = components_.result_converter->convert(
                clusters_for_read, read.read_id, read.len_raw_signal);
            for (const auto& result : results) {
                config_.result_writer->write(result);
            }
        }

        // Debug output to stdout
        output_ << read.read_id
                << "\tseeds=" << seeds_for_read.size()
                << "\thits=" << hits_for_read.size()
                << "\tanchors=" << clusters_for_read.expanded_anchor_count
                << "\tchains=" << clusters_for_read.clusters.size()
                << "\tlen=" << read.len_raw_signal
                << "\n";

        // Output all chains (one per line)
        for (std::size_t chain_idx = 0; chain_idx < clusters_for_read.clusters.size(); ++chain_idx) {
            const auto& chain = clusters_for_read.clusters[chain_idx];
            if (chain.anchors.empty()) continue;

            const auto& first_anchor = chain.anchors.front();
            const auto& last_anchor = chain.anchors.back();

            // Get path name
            std::string path_name;
            if (adj_store && first_anchor.path_id < adj_store->graph().pathCount()) {
                path_name = adj_store->graph().paths()[first_anchor.path_id].name;
            } else {
                path_name = "path" + std::to_string(first_anchor.path_id);
            }

            // Compute intervals (start of first anchor to end of last anchor)
            const auto q_start = first_anchor.read_pos;
            const auto q_end = last_anchor.read_pos + last_anchor.target.length;
            const auto r_start = first_anchor.ref_coord;
            const auto r_end = last_anchor.ref_coord + last_anchor.target.length;

            output_ << "  chain" << chain_idx << "\t"
                    << path_name << "\t"
                    << "score=" << static_cast<int>(chain.cluster_score) << "\t"
                    << "n=" << chain.anchors.size() << "\t"
                    << "q[" << q_start << "-" << q_end << "]\t"
                    << "r[" << r_start << "-" << r_end << "]\n";
        }

        // Blank line between reads
        output_ << "\n";
    }
}

}  // namespace piru::mapping
