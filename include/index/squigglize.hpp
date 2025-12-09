// SPDX-License-Identifier: MIT
// Squigglization and quantization helpers for indexing.

#pragma once

#include <vector>

#include "index/aln_graph.hpp"
#include "io/models/model.hpp"
#include "signal/alignment_quantizers/alignment_quantizer.hpp"
#include "signal/fuzzy_quantizers/fuzzy_quantizer.hpp"
#include "signal/signal_types.hpp"

namespace piru::index {

struct SquiggleResult {
    std::vector<std::vector<float>> raw_signals;
    std::vector<piru::signal::FuzzyQuantizedSignal> fuzzy_signals;
    std::vector<piru::signal::AlignmentQuantizedSignal> alignment_signals;
};

// Generate expected signals from the pore model and apply fuzzy/alignment quantization.
// - For each node: slide a k-mer window across the sequence and look up mean signal.
// - Normalize all samples globally (zero-mean, unit-std) before quantization.
// - Returns parallel vectors indexed by node ID.
SquiggleResult squigglizeAndQuantize(const AlnGraph& graph,
                                     const piru::io::KmerModel& model,
                                     const piru::signal::FuzzyQuantizer& fuzzy_quantizer,
                                     const piru::signal::AlignmentQuantizer& alignment_quantizer);

}  // namespace piru::index
