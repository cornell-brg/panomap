// SPDX-License-Identifier: MIT
// Interfaces for chaining backends.

#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "cli/parse.hpp"
#include "index/seed_store.hpp"

namespace piru::index {
class GraphStore;
}

namespace piru::mapping {

// Anchor in linear path coordinate space.
// Produced by expanding NodeAnchors using linearization coordinates.
// One NodeAnchor may expand to multiple PathAnchors (one per path occurrence).
struct PathAnchor {
    std::size_t query_pos{0};   // Position in query/read
    std::int64_t ref_coord{0};  // Linear position along reference path
    std::size_t length{0};      // Coverage length (from seed span)
    std::size_t path_id{0};     // Which reference path this anchor belongs to

    // Back-references for debugging and alignment construction
    std::size_t node_id{0};      // Graph node this anchor came from
    std::size_t node_offset{0};  // Offset within the node
};

// Minimal hit record used for clustering/chaining.
struct NodeAnchor {
    index::SeedEntry target;    // node_id + offset in graph
    std::size_t read_pos{0};  // seed position in the read
    std::uint64_t hash{0};    // seed hash (for debugging/uniqueness)
    std::size_t span{0};      // coverage length on query (from Seed.length, may be merged)
    std::size_t frequency{0};   // occurrences of this hash in the index
    mutable double score{0.0};  // Computed during clustering (mutable for legacy compatibility)
};

// Anchor candidate produced by clustering/chaining.
struct ChainedAnchor {
    index::SeedEntry target;      // node_id + offset in graph
    std::size_t read_pos{0};    // position in read
    double score{0.0};          // backend-specific score
    std::size_t chain_id{0};  // which chain this anchor belongs to

    // Optional: linear coordinates (for path-walk pipeline debugging)
    std::size_t path_id{0};     // reference path ID
    std::int64_t ref_coord{0};  // linear position on reference path
};

// A single chain: scored group of anchors.
struct Chain {
    double score{0.0};
    std::vector<ChainedAnchor> anchors;
};

struct ChainResult {
    double score{0.0};
    std::vector<ChainedAnchor> anchors;       // flat list from best chain
    std::vector<Chain> chains;             // all extracted chains
    std::size_t expanded_anchor_count{0};  // total anchors before chaining
};

// Abstract interface for chaining backends.
// Chainers select an optimal subset of seeds/anchors for alignment extension.
// The mapper treats this as a blackbox -- all backend-specific config and
// logic stays internal.
class Chainer {
public:
    virtual ~Chainer() = default;

    virtual ChainResult chain(const std::vector<PathAnchor>& anchors) const = 0;
    virtual std::string name() const = 0;

    // Debug: dump per-path chain diagnostics. Default no-op.
    virtual void dump_path_chains(const char* /*filename*/, const std::string& /*read_id*/,
                                  std::size_t /*read_length*/,
                                  const std::vector<PathAnchor>& /*anchors*/,
                                  const index::GraphStore* /*graph_store*/) const {}

    // Debug: dump per-anchor detail with chain membership. Default no-op.
    virtual void dump_anchor_detail(const char* /*filename*/, const std::string& /*read_id*/,
                                    std::size_t /*read_length*/,
                                    const std::vector<PathAnchor>& /*anchors*/,
                                    const index::GraphStore* /*graph_store*/) const {}
};

using ChainerPtr = std::unique_ptr<Chainer>;

// Factory: create a chainer from backend name + raw CLI args.
// Each backend parses what it needs from the CLI args.
ChainerPtr make_chainer(const std::string& backend, const cli::Parsed& parsed);

}  // namespace piru::mapping
