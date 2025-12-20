// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <vector>

#include "concurrency/executor.hpp"
#include "io/reads/read_provider.hpp"
#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
#include "signal/event_pipelines/event_pipeline_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"
#include "index/seed_store.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "mapping/anchor.hpp"
#include "mapping/seed_clusterer.hpp"
#include "mapping/anchor_expander.hpp"
#include "mapping/map_result.hpp"
#include "mapping/result_formatter.hpp"
#include "alignment/chain_aligner.hpp"
#include "io/results/result_writer.hpp"

namespace piru::mapping {

class SeedLookup {
public:
    SeedLookup(const index::SeedStore* store,
               const index::GraphStore* graph_store,
               std::size_t freq_threshold)
        : store_(store), graph_store_(graph_store), freq_threshold_(freq_threshold) {}

    void lookup(const signal::SeedBuffer& seeds, std::vector<SeedHitRecord>& out_hits) const;

private:
    const index::SeedStore* store_{nullptr};  // non-owning
    const index::GraphStore* graph_store_{nullptr};  // non-owning
    std::size_t freq_threshold_{0};
};

struct BatchMapperConfig {
    std::size_t batch_capacity_reads{4000};
    std::size_t batch_capacity_bytes{512 * 1024 * 1024};  // Reserved for future use.
    int num_threads{-1};  // -1 = automatic.

    signal::EventPipelineConfig event_pipeline_config{};  // Unified event detection + normalization
    signal::FuzzyQuantizerConfig fuzzy_config{};
    signal::AlignmentQuantizerConfig alignment_config{};
    signal::SeedExtractorConfig seed_config{};
    SeedClustererConfig clusterer_config{};
    const index::SeedStore* seed_store{nullptr};  // non-owning pointer to loaded SeedStore
    const index::GraphStore* graph_store{nullptr};  // non-owning pointer to loaded GraphStore
    const index::SignalStore* signal_store{nullptr};  // non-owning pointer to loaded SignalStore (for alignment)

    // Linearization coordinates (needed for DP chaining)
    // Non-owning pointer to linearization coords (from in-memory indexing or future deserialization)
    const std::vector<std::vector<index::LinearCoordinate>>* linearization_coords{nullptr};

    // Result writer for output (non-owning, optional)
    io::ResultWriter* result_writer{nullptr};

    // Alignment configuration (optional signal-level alignment)
    bool enable_alignment{false};
    alignment::ChainAlignerConfig align_config{};

    // Debug dump directories (empty = disabled)
    std::string dump_anchors_dir{};  // Dump anchors per read
    std::string dump_chains_dir{};   // Dump chains per read

    // Result formatting configuration
    ResultFormatterConfig formatter_config{};
};

struct BatchMapperStats {
    std::size_t batches{0};
    std::size_t reads_processed{0};
    std::size_t reads_mapped{0};      // reads with at least one chain
    std::size_t reads_unmapped{0};    // reads with no chains
    std::size_t results_written{0};   // total alignment results written
    std::size_t primary_alignments{0};
    std::size_t secondary_alignments{0};
};

struct BatchBuffer {
    std::vector<io::RawRead> raw_reads;
    std::vector<signal::NormalizedSignal> normalized;
    std::vector<signal::FuzzyQuantizedSignal> fuzzy_quantized;
    std::vector<signal::AlignmentQuantizedSignal> alignment_quantized;
    std::vector<signal::SeedBuffer> seeds;
    std::vector<std::vector<SeedHitRecord>> seed_hits;
    std::vector<std::vector<Anchor>> anchors;  // Debug: anchors after expansion
    std::vector<ReadMapResult> map_results;    // Unified mapping results
    std::size_t num_reads{0};

    void resize(std::size_t capacity);
    void clear();
};

struct PipelineComponents {
    signal::EventPipelinePtr event_pipeline;
    signal::FuzzyQuantizerPtr fuzzy_quantizer;
    signal::AlignmentQuantizerPtr alignment_quantizer;
    signal::SeedExtractorPtr seed_extractor;
    const index::SeedStore* seed_store{nullptr};  // non-owning; loaded index
    const index::GraphStore* graph_store{nullptr};  // non-owning; loaded index
    const index::SignalStore* signal_store{nullptr};  // non-owning; for alignment
    SeedLookup lookup{nullptr, nullptr, 0};
    AnchorExpanderPtr expander;  // Expands SeedHits to Anchors
    SeedClustererPtr clusterer;
    std::unique_ptr<ResultFormatter> result_formatter;       // Formats map results to PAF/GAF
    std::unique_ptr<alignment::ChainAligner> chain_aligner;  // Optional signal-level alignment
};

class BatchMapper {
public:
    BatchMapper(io::ReadProvider& provider, BatchMapperConfig config, std::ostream& output);

    BatchMapperStats process_all();

private:
    void load_batch(BatchBuffer& batch);
    void process_batch(BatchBuffer& batch);
    BatchMapperStats output_batch(const BatchBuffer& batch) const;
    void process_read(BatchBuffer& batch, std::size_t index) const;
    void lookup_seed_hits(const signal::SeedBuffer& seeds,
                          std::vector<SeedHitRecord>& hits_out) const;
    PipelineComponents create_components() const;

    BatchMapperConfig config_;
    io::ReadProvider& provider_;
    std::unique_ptr<concurrency::Executor> executor_;
    PipelineComponents components_;
    std::ostream& output_;
};

}  // namespace piru::mapping
