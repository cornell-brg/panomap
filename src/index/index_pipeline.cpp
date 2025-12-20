// SPDX-License-Identifier: MIT

#include "index/index_pipeline.hpp"

#include <memory>
#include <set>
#include <stdexcept>

#include "index/linearizer_factory.hpp"
#include "index/pseudo_linearize.hpp"
#include "index/seed_builder.hpp"
#include "index/squigglize.hpp"
#include "index/transform_dbg.hpp"
#include "index/vg_transform_factory.hpp"
#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
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

    IndexPipelineResult result;
    const std::size_t pore_k = model.k();

    // -------------------------------------------------------------------------
    // Stage 1: Graph Transformation (ImportedGraph -> AlnGraph)
    // -------------------------------------------------------------------------

    AlnGraph aln_graph;

    if (config.graph_flavor == "dbg") {
        if (config.graph_k == 0) {
            throw std::runtime_error("DBG graph requires graph_k to be set");
        }
        if (config.graph_k < pore_k) {
            throw std::runtime_error("graph k=" + std::to_string(config.graph_k) +
                                     " < pore k=" + std::to_string(pore_k) + " (invalid)");
        }

        aln_graph = transformDbg(imported, config.graph_k, pore_k);
    } else if (config.graph_flavor == "vg") {
        // VG transformation using path-guided approach
        TransformConfig transform_config;
        transform_config.uncovered_strategy = "expand";

        auto vg_transform = makeVGTransform("path_guided", transform_config);
        aln_graph = vg_transform->apply(imported, 0, pore_k);

        auto stats = vg_transform->getStats();
        LOG_INFO("VG transform: " + std::to_string(stats.original_node_count) +
                 " original nodes → " + std::to_string(stats.transformed_node_count) +
                 " transformed nodes (" + std::to_string(stats.node_expansion_ratio) +
                 "x expansion)");
        LOG_INFO("VG coverage: " + std::to_string(stats.uncovered_node_count) + " uncovered nodes");
    } else {
        throw std::runtime_error("Unknown graph flavor: " + config.graph_flavor);
    }

    if (!aln_graph.validate()) {
        throw std::runtime_error("AlnGraph validation failed after transformation");
    }

#ifdef PIRU_DUMP_GRAPHS
    GfaExporter::dumpAlnGraph(aln_graph, "transformed_graph.gfa", AlnGraphDumpMode::Bases);
#endif

    LOG_INFO("[1/4] transformed: " + std::to_string(aln_graph.nodeCount()) +
             " directional nodes (originally " + std::to_string(imported.nodes.size()) +
             " bidirected nodes)");

    // -------------------------------------------------------------------------
    // Stage 2: Linearization
    // -------------------------------------------------------------------------

    auto linearizer = make_linearizer(config.linearizer);
    if (!linearizer) {
        throw std::runtime_error("Failed to create linearizer: " + config.linearizer);
    }

    result.linearization_coords = linearizer->linearize(aln_graph);

    LOG_INFO("[2/4] linearized with " + linearizer->name() + " backend");

    // For superbubble backend, also store in graph nodes for serialization compatibility
    if (config.linearizer == "superbubble") {
        // Extract chain_id and linear_position from linearization coords
        for (std::size_t i = 0; i < aln_graph.nodeCount(); ++i) {
            if (!result.linearization_coords[i].empty()) {
                const auto& coord = result.linearization_coords[i][0];
                aln_graph.mutableNode(i).chain_id = static_cast<std::int64_t>(coord.path_id);
                aln_graph.mutableNode(i).linear_position = coord.ref_coord;
            } else {
                aln_graph.mutableNode(i).chain_id = -1;
                aln_graph.mutableNode(i).linear_position = -1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Stage 3: Squigglization + Quantization
    // -------------------------------------------------------------------------

    signal::FuzzyQuantizerConfig fuzzy_cfg;
    fuzzy_cfg.backend = config.fuzzy_quantizer;
    fuzzy_cfg.fine_min = config.fuzzy_fine_min;
    fuzzy_cfg.fine_max = config.fuzzy_fine_max;
    fuzzy_cfg.fine_range = config.fuzzy_fine_range;
    auto fuzzy_quantizer = signal::make_fuzzy_quantizer(fuzzy_cfg);
    if (!fuzzy_quantizer) {
        throw std::runtime_error("Failed to create fuzzy quantizer: " + config.fuzzy_quantizer);
    }

    signal::AlignmentQuantizerConfig align_cfg;
    align_cfg.backend = config.alignment_quantizer;
    if (config.alignment_scale > 0.0) {
        align_cfg.scale = config.alignment_scale;
    }
    auto alignment_quantizer = signal::make_alignment_quantizer(align_cfg);
    if (!alignment_quantizer) {
        throw std::runtime_error("Failed to create alignment quantizer: " +
                                 config.alignment_quantizer);
    }

    const auto squiggle_result = squigglizeAndQuantize(
        aln_graph, model, *fuzzy_quantizer, *alignment_quantizer);

#ifdef PIRU_DUMP_GRAPHS
    GfaExporter::dumpAlnGraph(aln_graph, "raw_signals.gfa", AlnGraphDumpMode::RawSignal,
                              &squiggle_result.raw_signals);
    GfaExporter::dumpAlnGraph(aln_graph, "fuzzy_quantized.gfa", AlnGraphDumpMode::FuzzyQuantized,
                              &squiggle_result.fuzzy_signals);
    GfaExporter::dumpAlnGraph(aln_graph, "aln_quantized.gfa", AlnGraphDumpMode::AlnQuantized,
                              &squiggle_result.alignment_signals);
#endif

    std::size_t total_samples = 0;
    std::set<std::int16_t> unique_tokens;
    for (const auto& sig : squiggle_result.fuzzy_signals) {
        total_samples += sig.tokens.size();
        for (const auto token : sig.tokens) {
            unique_tokens.insert(token);
        }
    }

    LOG_INFO("[3/4] squigglized: " + std::to_string(total_samples) + " signal samples, " +
             std::to_string(unique_tokens.size()) + " unique fuzzy tokens");

    // -------------------------------------------------------------------------
    // Stage 4: Seed Extraction & Indexing
    // -------------------------------------------------------------------------

    signal::SeedExtractorConfig extractor_cfg;
    extractor_cfg.backend = "kmer";
    extractor_cfg.k = config.seed_k;
    extractor_cfg.stride = config.seed_stride;
    extractor_cfg.qbits = 4;

    LOG_INFO("Seed extraction config: k=" + std::to_string(config.seed_k) +
             ", stride=" + std::to_string(config.seed_stride) +
             ", qbits=" + std::to_string(extractor_cfg.qbits) +
             ", backend=" + extractor_cfg.backend);

    auto extractor = signal::make_seed_extractor(extractor_cfg);
    if (!extractor) {
        throw std::runtime_error("Failed to create seed extractor");
    }

    SeedBuildConfig seed_cfg;
    seed_cfg.keep_least_frequent_fraction = config.seed_filter;

    auto seed_store = buildSeedStore(squiggle_result.fuzzy_signals, *extractor, seed_cfg);

    LOG_INFO("[4/4] indexed: " + std::to_string(seed_store.size()) +
             " unique seeds (max_freq=" + std::to_string(seed_store.max_hash_frequency()) + ")");

    // -------------------------------------------------------------------------
    // Package results
    // -------------------------------------------------------------------------

    result.graph_store = std::make_unique<AdjListGraphStore>(std::move(aln_graph));
    result.signal_store = std::make_unique<VectorSignalStore>(
        std::move(squiggle_result.alignment_signals));
    result.seed_store = std::make_unique<HashSeedStore>(std::move(seed_store));

    result.graph_flavor = imported.flavor;
    result.graph_k = config.graph_k;
    result.pore_k = pore_k;
    result.model_name = model.name();
    result.fuzzy_quantizer = fuzzy_cfg.backend;
    result.alignment_quantizer = align_cfg.backend;
    result.alignment_scale = alignment_quantizer->scale();
    result.alignment_offset = alignment_quantizer->offset();

    return result;
}

}  // namespace piru::index
