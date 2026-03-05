// SPDX-License-Identifier: MIT
// Shared indexing pipeline for building in-memory or serialized indexes.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/aln_graph.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/graphs/graph.hpp"
#include "io/models/model.hpp"

namespace piru::index {

// Configuration for the indexing pipeline.
//
// Single source of truth for indexing parameters (used by `piru index`).
struct IndexPipelineConfig {
    // -------------------------------------------------------------------------
    // Indexer Backend
    // -------------------------------------------------------------------------

    // Indexer backend: "node-first" or "path-walk"
    // - node-first: Two-pass approach processing node interiors once (faster for shared nodes)
    //               Uses global normalization across all nodes
    // - path-walk: Original per-path processing (per-path normalization)
    std::string indexer_backend{"path-walk"};

    // -------------------------------------------------------------------------
    // Signal Processing Parameters
    // -------------------------------------------------------------------------

    // Fuzzy quantizer: "rh2" (rawhash2) or "passthrough"
    // - Converts normalized signal to discrete tokens (default: 4-bit = 16 values)
    // - rh2: Adaptive quantization with fine/coarse regions
    std::string fuzzy_quantizer{"rh2"};
    float fuzzy_fine_min{-2.0f};    // Minimum value for fine quantization region
    float fuzzy_fine_max{2.0f};     // Maximum value for fine quantization region
    float fuzzy_fine_range{0.4f};   // Range per fine quantization bin
    std::uint32_t fuzzy_n_bins{0};  // Number of bins (0 = use 2^qbits = 16)

    // -------------------------------------------------------------------------
    // Seed Extraction Parameters
    // -------------------------------------------------------------------------

    // Seed extractor type: "kmer" or "minimizer"
    std::string seed_type{"minimizer"};

    // Minimizer window size (only used when seed_type="minimizer")
    // Within each window of w consecutive k-mers, only the minimum hash is kept.
    // w=1 recovers the rolling k-mer behavior.
    std::size_t minimizer_window{5};

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
    // - Default: 0.9 (configurable via CLI: --seed-freq-cutoff)
    double seed_freq_cutoff{0.9};

    // Subsample cap percentile for seeds above the frequency cutoff.
    // - When seed_freq_cutoff < 1.0, seeds above the threshold are subsampled
    //   down to the frequency at this percentile (instead of hard-dropped).
    // - Range: 0.0 to 1.0
    // - Default: 0.25 (configurable via CLI: --seed-freq-cap)
    double seed_freq_cap{0.25};

    // -------------------------------------------------------------------------
    // Debug Options
    // -------------------------------------------------------------------------

    // Dump per-path normalization stats to file (path-walk backend only)
    // Format: TSV with columns: path_name, mean, stddev, num_kmers
    std::string dump_norm_stats_path;

    // -------------------------------------------------------------------------
    // Parallelization
    // -------------------------------------------------------------------------

    // Executor for parallel indexing (optional).
    // If nullptr, indexing runs sequentially.
    // Caller owns the executor lifetime - must outlive the indexing call.
    concurrency::Executor* executor{nullptr};

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
    std::unique_ptr<SeedStore> seed_store;

    // Linearization coordinates (needed for DP chaining)
    std::vector<std::vector<LinearCoordinate>> linearization_coords;

    // Path lengths in base space (for anchor bounds checking)
    std::vector<std::size_t> path_lengths;

    // Metadata
    std::size_t pore_k{0};
    std::string model_name;
    std::string fuzzy_quantizer;
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
IndexPipelineResult run_index_pipeline(const io::ImportedGraph& imported,
                                       const io::KmerModel& model,
                                       const IndexPipelineConfig& config);

}  // namespace piru::index
