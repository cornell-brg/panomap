// SPDX-License-Identifier: MIT

#include "mapping/batch_mapper.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

// Dump hit statistics for a read to analyze frequency distributions.
void dumpHitStatsToFile(const char* filename,
                        const std::string& read_id,
                        const signal::SeedBuffer& seeds,
                        const std::vector<SeedHitRecord>& hits,
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
    out << "#SUMMARY\t" << read_id << "\t" << seeds.seeds.size() << "\t"
        << hashes_with_hits.size() << "\t" << total_hits << "\t" << freq_threshold << "\n";

    out << "#FREQ_STATS\tmin\tmax\tmedian\tp90\n";
    out << "#FREQ_STATS\t" << min_freq << "\t" << max_freq << "\t" << median_freq << "\t" << p90_freq << "\n";

    // Write per-seed details header
    out << "#SEEDS\thash\tread_pos\thit_count\n";

    // Group hits by hash to get per-seed hit counts
    std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> hash_info;  // hash -> (read_pos, freq)
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

// Dump best chain per path for diagnostic analysis.
// Runs DP chaining on anchors and extracts best chain per path by finding
// max dp[i] for each path_id and backtracking.
void dumpPathChainsToFile(const char* filename,
                          const std::string& read_id,
                          std::size_t read_length,
                          const std::vector<Anchor>& anchors,
                          const index::GraphStore* graph_store,
                          const SeedClustererConfig& clusterer_config) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        LOG_WARN("Failed to open path chains file: " + std::string(filename));
        return;
    }

    // Get path names if available
    const auto* adj_store = dynamic_cast<const index::AdjListGraphStore*>(graph_store);
    const std::vector<index::AlnPath>* paths = nullptr;
    std::size_t total_paths = 0;
    if (adj_store) {
        paths = &adj_store->graph().paths();
        total_paths = paths->size();
    }

    // Write header
    out << "read_id\tpath_id\tpath_name\tchain_score\tnum_anchors\t"
        << "query_start\tquery_end\tquery_span\tref_start\tref_end\tref_span\n";

    if (anchors.empty()) {
        // No anchors - output zeros for all paths
        for (std::size_t pid = 0; pid < total_paths; ++pid) {
            std::string path_name = paths ? (*paths)[pid].name : "path" + std::to_string(pid);
            out << read_id << "\t" << pid << "\t" << path_name << "\t"
                << "0\t0\t0\t0\t0.0\t0\t0\t0\n";
        }
        out.close();
        return;
    }

    // Sort anchors by (path_id, ref_coord, query_pos) - same as DPChainClusterer
    std::vector<Anchor> sorted_anchors = anchors;
    std::sort(sorted_anchors.begin(), sorted_anchors.end(),
              [](const Anchor& a, const Anchor& b) {
                  if (a.path_id != b.path_id) return a.path_id < b.path_id;
                  if (a.ref_coord != b.ref_coord) return a.ref_coord < b.ref_coord;
                  return a.query_pos < b.query_pos;
              });

    const std::size_t n = sorted_anchors.size();

    // DP arrays
    std::vector<double> dp(n, 0.0);
    std::vector<int> pred(n, -1);

    // DP parameters from config
    const double anchor_weight = clusterer_config.dp_anchor_weight;
    const double gap_penalty = clusterer_config.dp_gap_penalty;
    const double diag_penalty = clusterer_config.dp_diag_penalty;
    const double overlap_penalty = clusterer_config.dp_overlap_penalty;
    const std::size_t max_dist = clusterer_config.dp_max_dist;
    const std::size_t max_diag_dev = clusterer_config.dp_max_diag_dev;
    const std::size_t max_skip = clusterer_config.dp_max_skip;

    // Lambda: anchor score
    auto anchor_score = [anchor_weight](const Anchor& a) {
        return static_cast<double>(a.length) * anchor_weight;
    };

    // Lambda: can chain j -> i (query must strictly advance)
    auto can_chain = [max_dist, max_diag_dev](const Anchor& j, const Anchor& i) {
        if (i.path_id != j.path_id) return false;
        if (i.query_pos <= j.query_pos) return false;
        std::int64_t delta_ref = i.ref_coord - j.ref_coord;
        std::int64_t delta_query = static_cast<std::int64_t>(i.query_pos) -
                                   static_cast<std::int64_t>(j.query_pos);
        if (delta_ref < 0 || delta_query < 0) return false;
        if (static_cast<std::size_t>(delta_ref) > max_dist ||
            static_cast<std::size_t>(delta_query) > max_dist) return false;
        std::int64_t diag_dev = std::abs(delta_ref - delta_query);
        if (static_cast<std::size_t>(diag_dev) > max_diag_dev) return false;
        return true;
    };

    // Lambda: gap cost
    auto gap_cost = [gap_penalty, diag_penalty, overlap_penalty](const Anchor& j, const Anchor& i) {
        std::int64_t j_ref_end = j.ref_coord + static_cast<std::int64_t>(j.length);
        std::int64_t j_query_end = static_cast<std::int64_t>(j.query_pos + j.length);
        std::int64_t ref_gap = i.ref_coord - j_ref_end;
        std::int64_t query_gap = static_cast<std::int64_t>(i.query_pos) - j_query_end;
        double avg_gap = (std::max<std::int64_t>(0, ref_gap) + std::max<std::int64_t>(0, query_gap)) / 2.0;
        double cost = avg_gap * gap_penalty;
        std::int64_t delta_ref = i.ref_coord - j.ref_coord;
        std::int64_t delta_query = static_cast<std::int64_t>(i.query_pos) - static_cast<std::int64_t>(j.query_pos);
        cost += std::abs(delta_ref - delta_query) * diag_penalty;
        if (ref_gap < 0 || query_gap < 0) {
            double avg_overlap = (std::abs(std::min<std::int64_t>(0, ref_gap)) +
                                  std::abs(std::min<std::int64_t>(0, query_gap))) / 2.0;
            cost += avg_overlap * overlap_penalty;
        }
        return cost;
    };

    // DP loop - compute best score for each anchor
    for (std::size_t i = 0; i < n; ++i) {
        const auto& anchor_i = sorted_anchors[i];
        double best_score = anchor_score(anchor_i);
        int best_pred = -1;

        std::size_t num_skipped = 0;
        for (std::size_t j = i; j > 0 && num_skipped < max_skip; ) {
            --j;
            const auto& anchor_j = sorted_anchors[j];
            if (anchor_j.path_id != anchor_i.path_id) break;
            if (anchor_i.ref_coord - anchor_j.ref_coord > static_cast<std::int64_t>(max_dist)) break;
            if (!can_chain(anchor_j, anchor_i)) {
                ++num_skipped;
                continue;
            }
            num_skipped = 0;
            double cost = gap_cost(anchor_j, anchor_i);
            double score = dp[j] + anchor_score(anchor_i) - cost;
            if (score > best_score) {
                best_score = score;
                best_pred = static_cast<int>(j);
            }
        }
        dp[i] = best_score;
        pred[i] = best_pred;
    }

    // For each path, find anchor with max dp score and backtrack
    struct PathChainInfo {
        double score{0.0};
        std::vector<std::size_t> chain_indices;
    };
    std::unordered_map<std::size_t, PathChainInfo> best_per_path;

    // Find max dp score per path
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t pid = sorted_anchors[i].path_id;
        auto& info = best_per_path[pid];
        if (dp[i] > info.score) {
            info.score = dp[i];
            // Backtrack to get chain
            info.chain_indices.clear();
            int cur = static_cast<int>(i);
            while (cur != -1) {
                info.chain_indices.push_back(static_cast<std::size_t>(cur));
                cur = pred[cur];
            }
            std::reverse(info.chain_indices.begin(), info.chain_indices.end());
        }
    }

    // Collect all paths with their scores
    std::vector<std::pair<double, std::size_t>> scored_paths;
    for (std::size_t pid = 0; pid < total_paths; ++pid) {
        auto it = best_per_path.find(pid);
        double path_score = (it != best_per_path.end()) ? it->second.score : 0.0;
        scored_paths.emplace_back(path_score, pid);
    }
    // Sort by score descending
    std::sort(scored_paths.begin(), scored_paths.end(), std::greater<>());

    // Write one row per path
    for (const auto& [score, path_id] : scored_paths) {
        std::string path_name = "path" + std::to_string(path_id);
        if (paths && path_id < paths->size()) {
            path_name = (*paths)[path_id].name;
        }

        auto it = best_per_path.find(path_id);
        if (it == best_per_path.end() || it->second.chain_indices.empty()) {
            out << read_id << "\t" << path_id << "\t" << path_name << "\t"
                << "0\t0\t0\t0\t0.0\t0\t0\t0\n";
            continue;
        }

        const auto& chain_indices = it->second.chain_indices;

        // Compute query/ref spans from chain anchors
        std::size_t query_min = sorted_anchors[chain_indices.front()].query_pos;
        std::size_t query_max = query_min;
        std::int64_t ref_min = sorted_anchors[chain_indices.front()].ref_coord;
        std::int64_t ref_max = ref_min;

        for (std::size_t idx : chain_indices) {
            const auto& a = sorted_anchors[idx];
            query_min = std::min(query_min, a.query_pos);
            query_max = std::max(query_max, a.query_pos + a.length);
            ref_min = std::min(ref_min, a.ref_coord);
            ref_max = std::max(ref_max, a.ref_coord + static_cast<std::int64_t>(a.length));
        }

        std::size_t query_span = query_max - query_min;

        out << read_id << "\t"
            << path_id << "\t"
            << path_name << "\t"
            << std::fixed << std::setprecision(2) << it->second.score << "\t"
            << chain_indices.size() << "\t"
            << query_min << "\t"
            << query_max << "\t"
            << query_span << "\t"
            << ref_min << "\t"
            << ref_max << "\t"
            << (ref_max - ref_min) << "\n";
    }

    out.close();
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
        // Path-walk linearization → PathWalkExpander
        if (!config_.path_lengths) {
            throw std::runtime_error("PathWalkExpander requires path_lengths for bounds checking");
        }
        comps.expander = std::make_unique<PathWalkExpander>(
            *config_.linearization_coords, *config_.path_lengths);
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
            adj_store->graph(), nullptr, config_.formatter_config);
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

    // Dump hit stats if --dump-hit-stats is set
    static std::atomic<std::size_t> hit_stats_dump_counter{0};
    if (!config_.dump_hit_stats_dir.empty()) {
        std::size_t dump_num = hit_stats_dump_counter.fetch_add(1);
        std::string hit_stats_file = config_.dump_hit_stats_dir + "/read_" + std::to_string(dump_num) + "_hits.tsv";
        std::size_t freq_threshold = config_.seed_store ? config_.seed_store->frequency_threshold() : 0;
        dumpHitStatsToFile(hit_stats_file.c_str(), read.read_id, batch.seeds[index],
                           batch.seed_hits[index], freq_threshold);
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
        std::string anchor_file = config_.dump_anchors_dir + "/read_" + std::to_string(dump_num) + "_anchors.tsv";
        dumpAnchorsToFile(anchor_file.c_str(), anchors, read.read_id, config_.graph_store, nullptr);
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
        result.mappings.push_back(std::move(mapping));
    }

    // Dump chain if --dump-chains is set
    static std::atomic<std::size_t> chain_dump_counter{0};
    if (!config_.dump_chains_dir.empty()) {
        std::size_t dump_num = chain_dump_counter.fetch_add(1);
        std::string chain_file = config_.dump_chains_dir + "/read_" + std::to_string(dump_num) + "_chain.tsv";
        dumpChainToFile(chain_file.c_str(), result, read.read_id, config_.graph_store, nullptr);
    }

    // Dump best chain per path if --dump-path-chains is set (diagnostic)
    // Uses raw anchors (before clustering) to compute best chain per path via DP
    static std::atomic<std::size_t> path_chains_dump_counter{0};
    if (!config_.dump_path_chains_dir.empty()) {
        std::size_t dump_num = path_chains_dump_counter.fetch_add(1);
        std::string path_chains_file = config_.dump_path_chains_dir + "/read_" + std::to_string(dump_num) + "_path_chains.tsv";
        dumpPathChainsToFile(path_chains_file.c_str(), read.read_id, read.raw_signal.size(),
                             anchors, config_.graph_store, config_.clusterer_config);
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
