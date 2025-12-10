// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <vector>

#include "concurrency/executor.hpp"
#include "io/reads/read_provider.hpp"
#include "signal/alignment_quantizers/alignment_quantizer_factory.hpp"
#include "signal/event_detectors/event_detector_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/normalizers/normalizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"
#include "index/seed_store.hpp"

namespace piru::mapping {

struct SeedHitRecord {
    index::SeedHit target;      // node_id + offset in the graph
    std::size_t read_pos{0};    // position of the seed in the read
    std::uint64_t hash{0};      // seed hash to help downstream clustering/debugging
};

class SeedLookup {
public:
    SeedLookup(const index::SeedStore* store, std::size_t freq_threshold)
        : store_(store), freq_threshold_(freq_threshold) {}

    void lookup(const signal::SeedBuffer& seeds, std::vector<SeedHitRecord>& out_hits) const;

private:
    const index::SeedStore* store_{nullptr};  // non-owning
    std::size_t freq_threshold_{0};
};

struct BatchMapperConfig {
    std::size_t batch_capacity_reads{4000};
    std::size_t batch_capacity_bytes{512 * 1024 * 1024};  // Reserved for future use.
    int num_threads{-1};  // -1 = automatic.

    signal::EventDetectorConfig event_config{};
    signal::SignalNormalizerConfig normalizer_config{};
    signal::FuzzyQuantizerConfig fuzzy_config{};
    signal::AlignmentQuantizerConfig alignment_config{};
    signal::SeedExtractorConfig seed_config{};
    const index::SeedStore* seed_store{nullptr};  // non-owning pointer to loaded SeedStore
};

struct BatchMapperStats {
    std::size_t batches{0};
    std::size_t reads_processed{0};
};

struct BatchBuffer {
    std::vector<io::RawRead> raw_reads;
    std::vector<signal::EventSeries> events;
    std::vector<signal::NormalizedSignal> normalized;
    std::vector<signal::FuzzyQuantizedSignal> fuzzy_quantized;
    std::vector<signal::AlignmentQuantizedSignal> alignment_quantized;
    std::vector<signal::SeedBuffer> seeds;
    std::vector<std::vector<SeedHitRecord>> seed_hits;
    std::size_t num_reads{0};

    void resize(std::size_t capacity);
    void clear();
};

struct PipelineComponents {
    signal::EventDetectorPtr event_detector;
    signal::SignalNormalizerPtr normalizer;
    signal::FuzzyQuantizerPtr fuzzy_quantizer;
    signal::AlignmentQuantizerPtr alignment_quantizer;
    signal::SeedExtractorPtr seed_extractor;
    const index::SeedStore* seed_store{nullptr};  // non-owning; loaded index
    SeedLookup lookup{nullptr, 0};
};

class BatchMapper {
public:
    BatchMapper(io::ReadProvider& provider, BatchMapperConfig config, std::ostream& output);

    BatchMapperStats process_all();

private:
    void load_batch(BatchBuffer& batch);
    void process_batch(BatchBuffer& batch);
    void output_batch(const BatchBuffer& batch) const;
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
