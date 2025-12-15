// SPDX-License-Identifier: MIT
// Factory for seed clustering/chaining backends.

#include "mapping/seed_clusterer.hpp"

#include <cstdint>
#include <limits>
#include <cmath>
#include <unordered_map>

#include "mapping/dp_chain_clusterer.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Helper: select top seed per read position (legacy logic)
std::vector<const SeedHitRecord*> select_top_seeds_per_position(
    const std::vector<const SeedHitRecord*>& seeds) {
    if (seeds.empty()) return {};

    // Sort by read_pos ascending, with score/frequency/hash tie-breaking
    std::vector<const SeedHitRecord*> sorted = seeds;
    std::sort(sorted.begin(), sorted.end(),
              [](const SeedHitRecord* a, const SeedHitRecord* b) {
                  if (a->read_pos == b->read_pos) {
                      // Higher score first on ties
                      if (a->score == b->score) {
                          // Lower frequency preferred (more unique)
                          if (a->frequency == b->frequency) {
                              return a->hash < b->hash;  // Determinism
                          }
                          return a->frequency < b->frequency;
                      }
                      return a->score > b->score;
                  }
                  return a->read_pos < b->read_pos;
              });

    std::vector<const SeedHitRecord*> top;
    top.reserve(sorted.size());

    std::size_t current_pos = std::numeric_limits<std::size_t>::max();
    for (const auto* s : sorted) {
        if (s->read_pos != current_pos) {
            top.push_back(s);
            current_pos = s->read_pos;
        }
    }
    return top;
}

// Helper: compute probe stride (legacy logic)
std::size_t compute_probe_stride(std::size_t query_length, std::size_t max_probes) {
    if (max_probes == 0) return query_length == 0 ? 1 : query_length;
    std::size_t probes = std::max<std::size_t>(1, max_probes);
    if (query_length == 0) return 1;
    return (query_length + probes - 1) / probes;  // Ceil division
}

// Helper: sample probe seeds with stride (legacy logic)
std::vector<const SeedHitRecord*> sample_probe_seeds(
    const std::vector<const SeedHitRecord*>& seeds,
    std::size_t query_length,
    std::size_t max_probes,
    std::size_t stride_override) {
    if (seeds.empty() || max_probes == 0) return {};

    const std::size_t stride = (stride_override > 0)
                                   ? stride_override
                                   : compute_probe_stride(query_length, max_probes);

    // Reduce to top seed per read position
    auto top_per_pos = select_top_seeds_per_position(seeds);
    if (top_per_pos.empty()) return {};

    // Sample at stride intervals
    std::vector<const SeedHitRecord*> sampled;
    sampled.reserve(std::min(top_per_pos.size(), max_probes));

    std::size_t next_target = 0;
    for (const auto* s : top_per_pos) {
        if (sampled.empty()) {
            sampled.push_back(s);
            next_target = s->read_pos + stride;
            if (sampled.size() >= max_probes) break;
            continue;
        }
        if (s->read_pos >= next_target) {
            sampled.push_back(s);
            next_target = s->read_pos + stride;
            if (sampled.size() >= max_probes) break;
        }
    }

    // Ensure last position is considered if undersampled
    if (!top_per_pos.empty() && sampled.size() < max_probes) {
        const auto* last_seed = top_per_pos.back();
        if (sampled.empty() || sampled.back()->read_pos != last_seed->read_pos) {
            sampled.push_back(last_seed);
        }
    }

    // Cap to max_probes
    if (sampled.size() > max_probes) {
        sampled.resize(max_probes);
    }

    return sampled;
}

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

struct ScoredCluster {
    std::vector<const SeedHitRecord*> hits;  // All seeds in cluster
    std::size_t coverage{0};
    SeedAnchor top_seed{};
    double top_seed_score{0.0};
};

// Shared clustering logic (used by both FSE and Probe)
std::vector<ScoredCluster> cluster_and_score_seeds(
    const std::vector<SeedHitRecord>& hits,
    const SeedClustererConfig& config) {

    if (hits.empty()) return {};

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
            if (diag > last_diag + config.diagonal_cutoff) {
                // Finalize previous cluster if meets min size
                if (current.size() >= config.min_cluster_size) {
                    all_clusters.push_back(current);
                }
                current.clear();
            }
            current.push_back(h);
            last_diag = diag;
        }
        // Add final cluster if meets min size
        if (!current.empty() && current.size() >= config.min_cluster_size) {
            all_clusters.push_back(std::move(current));
        }
    }

    if (all_clusters.empty()) return {};

    // 3. Score ALL seeds in each cluster and find top seed
    std::vector<ScoredCluster> scored_clusters;
    scored_clusters.reserve(all_clusters.size());

    for (const auto& cluster_hits : all_clusters) {
        const std::size_t coverage = calculate_cluster_coverage(cluster_hits);

        // Score ALL seeds in this cluster (legacy behavior for probe sampling)
        double best_score = -std::numeric_limits<double>::infinity();
        SeedAnchor best_anchor{};

        for (const auto* h : cluster_hits) {
            const std::size_t uniqueness =
                (config.max_hash_frequency > h->frequency)
                    ? (config.max_hash_frequency - h->frequency)
                    : 0;
            const double score = static_cast<double>(coverage + uniqueness);

            // Populate score in the hit record (for probe sampling helpers)
            h->score = score;

            if (score > best_score) {
                best_score = score;
                best_anchor = SeedAnchor{
                    .target = h->target,
                    .read_pos = h->read_pos,
                    .score = score,
                    .cluster_id = 0,  // Will be set later after sorting
                };
            }
        }

        scored_clusters.push_back(ScoredCluster{
            .hits = cluster_hits,
            .coverage = coverage,
            .top_seed = best_anchor,
            .top_seed_score = best_score,
        });
    }

    // 4. Rank clusters by top seed score
    std::sort(scored_clusters.begin(), scored_clusters.end(),
              [](const ScoredCluster& a, const ScoredCluster& b) {
                  return a.top_seed_score > b.top_seed_score;
              });

    // 5. Assign cluster IDs based on rank
    for (std::size_t i = 0; i < scored_clusters.size(); ++i) {
        scored_clusters[i].top_seed.cluster_id = i;
    }

    // 6. Cap to max_clusters if specified
    if (config.max_clusters > 0 && static_cast<int>(scored_clusters.size()) > config.max_clusters) {
        scored_clusters.resize(static_cast<std::size_t>(config.max_clusters));
    }

    return scored_clusters;
}

class FseClusterer : public SeedClusterer {
public:
    explicit FseClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override {
        ClusterSummary summary;

        // For FSE/superbubble: 1:1 expansion (hits have chain_id already)
        summary.expanded_anchor_count = hits.size();

        // Use shared clustering logic
        auto scored_clusters = cluster_and_score_seeds(hits, config_);
        if (scored_clusters.empty()) return summary;

        // Return top seeds from selected clusters (multi-mapping support)
        summary.anchors.reserve(scored_clusters.size());
        for (const auto& cluster : scored_clusters) {
            summary.anchors.push_back(cluster.top_seed);
        }

        // Overall score is the best cluster score
        summary.score = scored_clusters[0].top_seed_score;

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
        summary.expanded_anchor_count = hits.size();
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

class ProbeClusterer : public SeedClusterer {
public:
    explicit ProbeClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<SeedHitRecord>& hits) const override {
        ClusterSummary summary;

        // For Probe/superbubble: 1:1 expansion (hits have chain_id already)
        summary.expanded_anchor_count = hits.size();

        // Use shared clustering logic (same as FSE)
        auto scored_clusters = cluster_and_score_seeds(hits, config_);
        if (scored_clusters.empty()) return summary;

        // Estimate query length from seed hits (for probe stride calculation)
        std::size_t query_length = 0;
        for (const auto& h : hits) {
            query_length = std::max(query_length, h.read_pos + h.span);
        }
        if (query_length == 0) query_length = 1000;  // Fallback if no hits

        // Sample probes from each cluster and return as grouped clusters
        summary.clusters.reserve(scored_clusters.size());

        for (const auto& scored_cluster : scored_clusters) {
            // Sample probe seeds from this cluster
            auto probes = sample_probe_seeds(
                scored_cluster.hits,
                query_length,
                config_.max_probes_per_cluster,
                config_.probe_stride);

            if (probes.empty()) continue;

            // Build ClusterGroup with sampled probes
            ClusterGroup group;
            group.cluster_score = scored_cluster.top_seed_score;
            group.anchors.reserve(probes.size());

            for (const auto* probe : probes) {
                group.anchors.push_back(SeedAnchor{
                    .target = probe->target,
                    .read_pos = probe->read_pos,
                    .score = scored_cluster.top_seed_score,  // Use cluster score
                    .cluster_id = scored_cluster.top_seed.cluster_id,
                });
            }

            summary.clusters.push_back(std::move(group));
        }

        // Overall score is the best cluster score
        if (!scored_clusters.empty()) {
            summary.score = scored_clusters[0].top_seed_score;
        }

        // Also populate flat anchors list for compatibility (all probes from all clusters)
        for (const auto& group : summary.clusters) {
            summary.anchors.insert(summary.anchors.end(),
                                   group.anchors.begin(), group.anchors.end());
        }

        return summary;
    }

    std::string name() const override { return "probe"; }

private:
    SeedClustererConfig config_;
};

}  // namespace

SeedClustererPtr make_seed_clusterer(
    const SeedClustererConfig& config,
    const std::vector<std::vector<index::LinearCoordinate>>* linearization_coords,
    const index::GraphStore* graph_store) {

    if (config.backend == "fse" || config.backend.empty()) {
        return std::make_unique<FseClusterer>(config);
    }
    if (config.backend == "probe") {
        return std::make_unique<ProbeClusterer>(config);
    }
    if (config.backend == "dp-chain") {
        if (!linearization_coords) {
            LOG_ERROR("DP chaining requires linearization_coords (use --linearizer path-walk)");
            return std::make_unique<NoopClusterer>();
        }
        // Use default DPChainClustererConfig for MVP (parameters can be exposed via CLI later)
        DPChainClustererConfig dp_config;
        dp_config.min_chain_score = 0;  // MVP: accept any chain for testing
        return std::make_unique<DPChainClusterer>(*linearization_coords, dp_config);
    }
    if (config.backend == "chaining" || config.backend == "noop") {
        // Placeholder: chaining to be implemented; fallback to noop for now.
        LOG_WARN("Seed clusterer backend '" + config.backend + "' not implemented, using noop");
        return std::make_unique<NoopClusterer>();
    }
    LOG_WARN("Unknown seed clusterer backend '" + config.backend + "', using noop");
    return std::make_unique<NoopClusterer>();
}

}  // namespace piru::mapping
