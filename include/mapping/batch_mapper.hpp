// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "concurrency/executor.hpp"
#include "index/graph_store.hpp"
#include "index/linearizer.hpp"
#include "index/seed_store.hpp"
#include "io/reads/read_provider.hpp"
#include "io/results/result_writer.hpp"
#include "mapping/chainer.hpp"
#include "mapping/map_result.hpp"
#include "mapping/result_formatter.hpp"
#include "signal/event_pipelines/event_pipeline_factory.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer_factory.hpp"
#include "signal/seed_extractors/seed_extractor_factory.hpp"
#include "signal/signal_types.hpp"

namespace piru::mapping {

class SeedLookup {
public:
  SeedLookup(const index::SeedStore* store, std::size_t freq_threshold)
      : store_(store), freq_threshold_(freq_threshold) {}

  void lookup(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& out_hits) const;

private:
  const index::SeedStore* store_{nullptr};  // non-owning
  std::size_t freq_threshold_{0};
};

struct BatchMapperConfig {
  std::size_t batch_capacity_reads{4000};
  std::size_t batch_capacity_bytes{512 * 1024 * 1024};  // Reserved for future use.
  int num_threads{-1};                                  // -1 = automatic.

  signal::EventPipelineConfig event_pipeline_config{};  // Unified event detection + normalization
  signal::FuzzyQuantizerConfig fuzzy_config{};
  signal::SeedExtractorConfig seed_config{};
  std::string chainer_backend{"path-chain"};
  cli::Parsed chainer_parsed{};                   // CLI args forwarded to chainer
  const index::SeedStore* seed_store{nullptr};    // non-owning pointer to loaded SeedStore
  const index::GraphStore* graph_store{nullptr};  // non-owning pointer to loaded GraphStore
  // Linearization coordinates (needed for DP chaining)
  // Non-owning pointer to linearization coords (from in-memory indexing or future
  // deserialization)
  const std::vector<std::vector<index::LinearCoordinate>>* linearization_coords{nullptr};
  // Path lengths for anchor bounds checking (parallel to graph paths)
  const std::vector<std::size_t>* path_lengths{nullptr};

  // Result writer for output (non-owning, optional)
  io::ResultWriter* result_writer{nullptr};

  // Anchor merging (passed to chainer)
  bool enable_anchor_merge{true};

  // ROI classification
  const std::unordered_set<std::size_t>* roi_nodes{nullptr};  // non-owning
  std::string classify_mode{};  // "enrich" or "deplete", empty = disabled
  bool roi_filter_anchors{false};     // --chain-target: filter anchors to ROI nodes
  double roi_score_threshold{30.0};   // --chain-target: min chain score for accept
  double roi_overlap_threshold{0.5};  // --chain-genome: min ROI overlap for accept

  // Result formatting configuration
  ResultFormatterConfig formatter_config{};
};

struct BatchMapperStats {
  std::size_t batches{0};
  std::size_t reads_processed{0};
  std::size_t reads_mapped{0};     // reads with at least one chain
  std::size_t reads_unmapped{0};   // reads with no chains
  std::size_t results_written{0};  // total alignment results written
  std::size_t primary_alignments{0};
  std::size_t secondary_alignments{0};
};

struct BatchBuffer {
  std::vector<io::RawRead> raw_reads;
  std::vector<signal::NormalizedSignal> normalized;
  std::vector<signal::FuzzyQuantizedSignal> fuzzy_quantized;
  std::vector<signal::SeedBuffer> seeds;
  std::vector<std::vector<NodeAnchor>> seed_hits;
  std::vector<ReadMapResult> map_results;
  std::size_t num_reads{0};

  void resize(std::size_t capacity);
  void clear();
};

struct PipelineComponents {
  signal::EventPipelinePtr event_pipeline;
  signal::FuzzyQuantizerPtr fuzzy_quantizer;
  signal::SeedExtractorPtr seed_extractor;
  const index::SeedStore* seed_store{nullptr};    // non-owning; loaded index
  const index::GraphStore* graph_store{nullptr};  // non-owning; loaded index
  SeedLookup lookup{nullptr, 0};
  ChainerPtr chainer;
  std::unique_ptr<ResultFormatter> result_formatter;  // Formats map results to PAF/GAF
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
  void lookup_seed_hits(const signal::SeedBuffer& seeds, std::vector<NodeAnchor>& hits_out) const;
  PipelineComponents create_components() const;

  BatchMapperConfig config_;
  io::ReadProvider& provider_;
  std::unique_ptr<concurrency::Executor> executor_;
  PipelineComponents components_;
  std::ostream& output_;
};

}  // namespace piru::mapping
