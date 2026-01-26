// SPDX-License-Identifier: MIT

#include "index/index_pipeline.hpp"

#include <chrono>
#include <memory>
#include <stdexcept>

#include "index/node_first_indexer.hpp"
#include "index/path_walk_indexer.hpp"
#include "index/simple_expand.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "util/logging.hpp"

#ifdef PIRU_DUMP_GRAPHS
#include "io/graphs/gfa_exporter.hpp"
#endif

namespace piru::index {

IndexPipelineResult run_index_pipeline(
    const io::ImportedGraph& imported,
    const io::KmerModel& model,
    const IndexPipelineConfig& config) {

    auto stage_start = std::chrono::high_resolution_clock::now();

    // Stage 1: Simple ±expand (2x nodes)
    AlnGraph aln_graph = simpleExpand(imported);

    auto stage_elapsed = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - stage_start).count();
    LOG_INFO("[1/2] Transforming graph to directional graph: " + std::to_string(aln_graph.nodeCount()) +
             " nodes, " + std::to_string(aln_graph.edgeCount()) + " edges, " +
             std::to_string(aln_graph.pathCount()) + " paths [" +
             std::to_string(stage_elapsed) + "s]");

#ifdef PIRU_DUMP_GRAPHS
    GfaExporter::dumpAlnGraph(aln_graph, "simple_expanded.gfa", AlnGraphDumpMode::Bases);
    LOG_INFO("Exported simple expanded graph to simple_expanded.gfa");
#endif

    // Stage 2: Indexing (squigglize + linearize + seed extraction)
    stage_start = std::chrono::high_resolution_clock::now();

    // Create fuzzy quantizer
    signal::FuzzyQuantizerConfig fuzzy_cfg;
    fuzzy_cfg.backend = config.fuzzy_quantizer;
    fuzzy_cfg.pore_model = model.name();
    fuzzy_cfg.fine_min = config.fuzzy_fine_min;
    fuzzy_cfg.fine_max = config.fuzzy_fine_max;
    fuzzy_cfg.fine_range = config.fuzzy_fine_range;
    fuzzy_cfg.n_bins = config.fuzzy_n_bins;
    auto fuzzy_quantizer = signal::make_fuzzy_quantizer(fuzzy_cfg);
    if (!fuzzy_quantizer) {
        throw std::runtime_error("Failed to create fuzzy quantizer: " + config.fuzzy_quantizer);
    }

    // Create seed extractor
    signal::SeedExtractorConfig extractor_cfg;
    extractor_cfg.backend = "kmer";
    extractor_cfg.k = config.seed_k;
    extractor_cfg.stride = config.seed_stride;
    extractor_cfg.qbits = 4;
    auto extractor = signal::make_seed_extractor(extractor_cfg);
    if (!extractor) {
        throw std::runtime_error("Failed to create seed extractor");
    }

    // Package result - populated by either backend
    IndexPipelineResult result;
    std::vector<std::size_t> path_lengths;

    if (config.indexer_backend == "node-first") {
        // Node-first indexing: process node interiors once, then fill boundaries
        NodeFirstIndexConfig nfi_config;
        nfi_config.seed_k = config.seed_k;
        nfi_config.seed_stride = config.seed_stride;
        nfi_config.seed_filter = config.seed_filter;
        nfi_config.executor = config.executor;

        auto nfi_result = nodeFirstIndex(aln_graph, model, *fuzzy_quantizer, *extractor, nfi_config);

        stage_elapsed = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - stage_start).count();
        LOG_INFO("[2/2] node-first indexed: " + std::to_string(nfi_result.seeds_unique) +
                 " unique seeds (global_mean=" + std::to_string(nfi_result.global_mean) +
                 ", global_std=" + std::to_string(nfi_result.global_std) +
                 ") [" + std::to_string(stage_elapsed) + "s]");

        path_lengths = std::move(nfi_result.path_lengths);
        result.seed_store = std::move(nfi_result.seed_store);
        result.linearization_coords = std::move(nfi_result.linearization_coords);
    } else {
        // Path-walk indexing: process each path independently (per-path normalization)
        PathWalkIndexConfig pwi_config;
        pwi_config.seed_k = config.seed_k;
        pwi_config.seed_stride = config.seed_stride;
        pwi_config.seed_filter = config.seed_filter;
        pwi_config.dump_norm_stats_path = config.dump_norm_stats_path;
        pwi_config.executor = config.executor;

        auto pwi_result = pathWalkIndex(aln_graph, model, *fuzzy_quantizer, *extractor, pwi_config);

        stage_elapsed = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - stage_start).count();
        LOG_INFO("[2/2] path-walk indexed: " + std::to_string(pwi_result.seeds_unique) +
                 " unique seeds [" + std::to_string(stage_elapsed) + "s]");

        path_lengths = std::move(pwi_result.path_lengths);
        result.seed_store = std::move(pwi_result.seed_store);
        result.linearization_coords = std::move(pwi_result.linearization_coords);
    }

    // Copy path lengths to graph paths (for result_formatter coordinate flipping)
    for (std::size_t i = 0; i < aln_graph.pathCount(); ++i) {
        aln_graph.mutablePath(i).length = path_lengths[i];
    }

    // Package remaining result fields
    result.graph_store = std::make_unique<AdjListGraphStore>(std::move(aln_graph));
    result.path_lengths = std::move(path_lengths);  // For anchor bounds checking
    result.pore_k = model.k();
    result.model_name = model.name();
    result.fuzzy_quantizer = fuzzy_cfg.backend;

    return result;
}

}  // namespace piru::index
