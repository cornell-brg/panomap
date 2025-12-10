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

// Helper to calculate coverage with interval merging (legacy logic)
std::size_t calculate_cluster_coverage(const std::vector<const SeedHitRecord*>& cluster) {
    if (cluster.empty()) return 0;

    // Build coverage intervals [start, end]
    std::vector<std::pair<std::size_t, std::size_t>> intervals;
    intervals.reserve(cluster.size());
    for (const auto* h : cluster) {
        std::size_t start = h->read_pos;
        std::size_t end = h->read_pos + h->span - 1;
        intervals.emplace_back(start, end);
    }

    // Sort and merge overlapping intervals
    std::sort(intervals.begin(), intervals.end());

    std::size_t total_coverage = 0;
    std::size_t current_start = intervals[0].first;
    std::size_t current_end = intervals[0].second;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        auto [start, end] = intervals[i];
        if (start <= current_end + 1) {
            // Overlapping or adjacent — extend
            current_end = std::max(current_end, end);
        } else {
            // No overlap — finalize current interval
            total_coverage += current_end - current_start + 1;
            current_start = start;
            current_end = end;
        }
    }

    // Add the last interval
    total_coverage += current_end - current_start + 1;
    return total_coverage;
}

struct ClusterRep {
    std::size_t coverage{0};
    SeedAnchor top_seed{};
    double top_seed_score{0.0};
};

class FseClusterer : public SeedClusterer {
public:
    explicit FseClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override {
        ClusterSummary summary;
        if (hits.empty()) return summary;

        // 1. Group seeds by chain_id (legacy logic: ignore seeds without chain)
        std::unordered_map<std::int64_t, std::vector<const SeedHitRecord*>> by_chain;
        by_chain.reserve(hits.size());
        for (const auto& h : hits) {
            if (!h.chain_id.has_value()) continue;  // Skip non-chained seeds
            by_chain[*h.chain_id].push_back(&h);
        }

        // 2. Within each chain, cluster by diagonal with gap tolerance
        std::vector<std::vector<const SeedHitRecord*>> all_clusters;

        for (auto& [chain_id, chain_hits] : by_chain) {
            // Sort by diagonal (GA-style: ref_pos - read_pos)
            std::sort(chain_hits.begin(), chain_hits.end(),
                      [](const SeedHitRecord* a, const SeedHitRecord* b) {
                          const std::int64_t da = a->linear_pos.value() - static_cast<std::int64_t>(a->read_pos);
                          const std::int64_t db = b->linear_pos.value() - static_cast<std::int64_t>(b->read_pos);
                          if (da == db) return a->read_pos < b->read_pos;  // Tie-break
                          return da < db;
                      });

            // Split clusters when diagonal gap > cutoff
            std::vector<const SeedHitRecord*> current;
            current.reserve(chain_hits.size());
            std::int64_t last_diag = 0;
            bool first = true;

            for (const auto* h : chain_hits) {
                const std::int64_t diag = h->linear_pos.value() - static_cast<std::int64_t>(h->read_pos);
                if (first) {
                    current.push_back(h);
                    last_diag = diag;
                    first = false;
                    continue;
                }
                if (diag > last_diag + config_.diagonal_cutoff) {
                    // Finalize previous cluster if meets min size
                    if (current.size() >= config_.min_cluster_size) {
                        all_clusters.push_back(current);
                    }
                    current.clear();
                }
                current.push_back(h);
                last_diag = diag;
            }
            // Add final cluster if meets min size
            if (!current.empty() && current.size() >= config_.min_cluster_size) {
                all_clusters.push_back(std::move(current));
            }
        }

        if (all_clusters.empty()) return summary;

        // 3. Score seeds and find top seed per cluster
        std::vector<ClusterRep> cluster_reps;
        cluster_reps.reserve(all_clusters.size());

        for (const auto& cluster_hits : all_clusters) {
            const std::size_t coverage = calculate_cluster_coverage(cluster_hits);

            // Find best seed in this cluster
            double best_score = -std::numeric_limits<double>::infinity();
            SeedAnchor best_anchor{};

            for (const auto* h : cluster_hits) {
                const std::size_t uniqueness =
                    (config_.max_hash_frequency > h->frequency)
                        ? (config_.max_hash_frequency - h->frequency)
                        : 0;
                const double score = static_cast<double>(coverage + uniqueness);

                if (score > best_score) {
                    best_score = score;
                    best_anchor = SeedAnchor{
                        .target = h->target,
                        .read_pos = h->read_pos,
                        .score = score,
                    };
                }
            }

            cluster_reps.push_back(ClusterRep{
                .coverage = coverage,
                .top_seed = best_anchor,
                .top_seed_score = best_score,
            });
        }

        // 4. Rank clusters by top seed score
        std::sort(cluster_reps.begin(), cluster_reps.end(),
                  [](const ClusterRep& a, const ClusterRep& b) {
                      return a.top_seed_score > b.top_seed_score;
                  });

        // 5. Cap to max_clusters if specified
        if (config_.max_clusters > 0 && static_cast<int>(cluster_reps.size()) > config_.max_clusters) {
            cluster_reps.resize(static_cast<std::size_t>(config_.max_clusters));
        }

        // 6. Return top seeds from selected clusters (multi-mapping support)
        summary.anchors.reserve(cluster_reps.size());
        for (const auto& rep : cluster_reps) {
            summary.anchors.push_back(rep.top_seed);
        }

        // Overall score is the best cluster score
        if (!cluster_reps.empty()) {
            summary.score = cluster_reps[0].top_seed_score;
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
