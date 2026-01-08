// SPDX-License-Identifier: MIT
// Shared indexing pipeline for building in-memory or serialized indexes.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "index/aln_graph.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "index/signal_store.hpp"
#include "io/graphs/graph.hpp"
#include "io/models/model.hpp"

namespace piru::index {

// Configuration for the indexing pipeline.
//
// IMPORTANT: These defaults are shared by both `piru index` and `piru map --graph`.
// This is the single source of truth for indexing parameters.
struct IndexPipelineConfig {
    // -------------------------------------------------------------------------
    // Pipeline Mode
    // -------------------------------------------------------------------------

    // Pipeline backend: "simple" or "classic"
    // - simple: all in one go path indexing (current default)
    // - classic: Existing path-guided transform with context handling
    std::string pipeline_mode{"simple"};

    // -------------------------------------------------------------------------
    // Linearization Parameters
    // -------------------------------------------------------------------------

    // Linearization backend: "superbubble" or "path-walk"
    // - superbubble: Uses pseudo-linearization with local chain coordinates
    //                (supports serialization, works without reference paths)
    // - path-walk:   Walks reference paths to assign global coordinates
    //                (requires graph with paths, needed for DP chaining)
    std::string linearizer{"superbubble"};

    // -------------------------------------------------------------------------
    // Signal Processing Parameters
    // -------------------------------------------------------------------------

    // Fuzzy quantizer: "rh2" (rawhash2) or "passthrough"
    // - Converts normalized signal to discrete tokens (default: 4-bit = 16 values)
    // - rh2: Adaptive quantization with fine/coarse regions
    std::string fuzzy_quantizer{"rh2"};
    float fuzzy_fine_min{-2.0f};   // Minimum value for fine quantization region
    float fuzzy_fine_max{2.0f};    // Maximum value for fine quantization region
    float fuzzy_fine_range{0.4f};  // Range per fine quantization bin
    std::uint32_t fuzzy_n_bins{0}; // Number of bins (0 = use 2^qbits = 16)

    // Alignment quantizer: "int16", "int8", or "passthrough"
    // - Converts normalized signal to integer format for alignment
    // - int16: 16-bit signed integers with auto-scaling (default)
    std::string alignment_quantizer{"int16"};

    // Alignment quantizer scale override (0 = auto-detect from signal range)
    double alignment_scale{0.0};

    // -------------------------------------------------------------------------
    // Seed Extraction Parameters
    // -------------------------------------------------------------------------

    // Seed extraction mode: "node" or "path"
    // - node: Extract seeds from each node independently (current behavior)
    //         Seeds cannot cross node boundaries
    // - path: Walk paths and extract seeds that can cross node boundaries
    //         Better coverage for graphs with short nodes (e.g., VG-built)
    //         Deduplicates seeds from shared regions across paths
    std::string seed_mode{"path"};

    // Seed k-mer size (number of fuzzy tokens hashed together)
    // - Larger k = more specific seeds, fewer false matches
    // - Smaller k = more sensitive but noisier matches
    // - Default: 6 tokens per seed
    std::size_t seed_k{6};

    // Seed extraction stride (spacing between consecutive seeds)
    // - stride=1: Extract seed at every position (dense)
    // - stride>1: Sample every Nth position (sparse, faster)
    // - Default: 1 (extract all positions)
    std::size_t seed_stride{1};

    // Seed frequency filter (fraction of least frequent seeds to keep)
    // - Range: 0.0 to 1.0
    // - 0.9 = keep bottom 90% least frequent seeds (filter top 10% repetitive)
    // - 1.0 = keep all seeds (no filtering)
    // - Lower values reduce index size and mapping noise
    // - Default: 0.9 (configurable via CLI: --seed-filter)
    double seed_filter{0.5};

    // -------------------------------------------------------------------------
    // Note on Additional Parameters
    // -------------------------------------------------------------------------
    //
    // Fuzzy quantizer qbits (vocabulary size):
    //   - Hardcoded to 4 bits in fuzzy_quantizer.hpp (16 possible token values)
    //   - Seed extractor qbits MUST match (hardcoded to 4 in index_pipeline.cpp)
    //   - Do not change without coordinating both values
    //
    // Pore model k-mer size:
    //   - Determined by the selected model (e.g., r9.4 → k=6, r10.4 → k=9)
    //   - Not configurable (model-dependent)
};

// Result of running the indexing pipeline (in-memory representation).
struct IndexPipelineResult {
    // Core index components
    std::unique_ptr<GraphStore> graph_store;
    std::unique_ptr<SignalStore> signal_store;
    std::unique_ptr<SeedStore> seed_store;

    // Linearization coordinates (needed for DP chaining)
    std::vector<std::vector<LinearCoordinate>> linearization_coords;

    // Metadata
    std::size_t pore_k{0};
    std::string model_name;
    std::string fuzzy_quantizer;
    std::string alignment_quantizer;
    double alignment_scale{1.0};
    double alignment_offset{0.0};
};

// Run the full indexing pipeline on an imported graph.
//
// Pipeline stages:
// 1. Transform: ImportedGraph → AlnGraph (VG path-guided)
// 2. Linearize: Assign coordinates to nodes (path-walk)
// 3. Squigglize: Generate reference signals for each node
// 4. Quantize: Convert signals to fuzzy and alignment quantized forms
// 5. Build seeds: Extract and index k-mer seeds from fuzzy signals
//
// Returns an IndexPipelineResult with all components in memory.
IndexPipelineResult run_index_pipeline(
    const io::ImportedGraph& imported,
    const io::KmerModel& model,
    const IndexPipelineConfig& config);

}  // namespace piru::index
