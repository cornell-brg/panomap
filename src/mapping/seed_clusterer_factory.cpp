// SPDX-License-Identifier: MIT
// Factory for seed clustering/chaining backends.

#include "mapping/seed_clusterer.hpp"

#include <cstdint>
#include <limits>
#include <cmath>
#include <unordered_map>

#include "util/logging.hpp"

namespace piru::mapping {

namespace {

struct ClusterKey {
    std::int64_t chain_id;
    std::int64_t diag;
    bool operator==(const ClusterKey& other) const {
        return chain_id == other.chain_id && diag == other.diag;
    }
};

struct ClusterKeyHash {
    std::size_t operator()(const ClusterKey& k) const {
        // Simple combine
        return std::hash<std::int64_t>{}(k.chain_id) ^ (std::hash<std::int64_t>{}(k.diag) << 1);
    }
};

class FseClusterer : public SeedClusterer {
public:
    explicit FseClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override {
        ClusterSummary summary;
        if (hits.empty()) return summary;

        std::unordered_map<ClusterKey, std::vector<const SeedHitRecord*>, ClusterKeyHash> clusters;
        clusters.reserve(hits.size());

        for (const auto& h : hits) {
            const std::int64_t chain = h.chain_id.has_value() ? *h.chain_id : -1;
            const std::int64_t lin = h.linear_pos.has_value() ? *h.linear_pos : 0;
            const std::int64_t diag = lin - static_cast<std::int64_t>(h.read_pos);
            clusters[{chain, diag}].push_back(&h);
        }

        double best_score = -std::numeric_limits<double>::infinity();
        SeedAnchor best_anchor{};

        for (const auto& [key, cluster_hits] : clusters) {
            std::size_t coverage = 0;
            for (const auto* h : cluster_hits) {
                coverage += h->span;
            }
            for (const auto* h : cluster_hits) {
                const std::size_t uniq =
                    (config_.max_hash_frequency > h->frequency)
                        ? (config_.max_hash_frequency - h->frequency)
                        : 0;
                const double score = static_cast<double>(coverage + uniq);
                if (score > best_score) {
                    best_score = score;
                    best_anchor = SeedAnchor{
                        .target = h->target,
                        .read_pos = h->read_pos,
                        .score = score,
                    };
                }
            }
        }

        if (std::isfinite(best_score)) {
            summary.score = best_score;
            summary.anchors.push_back(best_anchor);
        }
        return summary;
    }

    std::string name() const override { return "fse"; }

private:
    SeedClustererConfig config_;
};

class NoopClusterer : public SeedClusterer {
public:
    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override {
        ClusterSummary summary;
        summary.anchors.reserve(hits.size());
        for (const auto& h : hits) {
            summary.anchors.push_back(SeedAnchor{
                .target = h.target,
                .read_pos = h.read_pos,
                .score = 0.0,
            });
        }
        return summary;
    }

    std::string name() const override { return "noop"; }
};

}  // namespace

SeedClustererPtr make_seed_clusterer(const SeedClustererConfig& config) {
    if (config.backend == "fse" || config.backend.empty()) {
        return std::make_unique<FseClusterer>(config);
    }
    if (config.backend == "probe" || config.backend == "chaining" || config.backend == "noop") {
        // Placeholder: probe/chaining to be implemented; fallback to noop for now.
        LOG_WARN("Seed clusterer backend '" + config.backend + "' not implemented, using noop");
        return std::make_unique<NoopClusterer>();
    }
    LOG_WARN("Unknown seed clusterer backend '" + config.backend + "', using noop");
    return std::make_unique<NoopClusterer>();
}

}  // namespace piru::mapping
