// SPDX-License-Identifier: MIT
// Factory for seed clustering/chaining backends.

#include "mapping/seed_clusterer.hpp"

#include <cstdint>
#include <limits>
#include <cmath>
#include <unordered_map>

#include "mapping/anchor.hpp"
#include "mapping/dp_chain_clusterer.hpp"
#include "util/logging.hpp"

namespace piru::mapping {

namespace {

// Helper: select top anchor per query position (for probe sampling)
std::vector<const Anchor*> select_top_anchors_per_position(
    const std::vector<const Anchor*>& anchors) {
    if (anchors.empty()) return {};

    // Sort by query_pos ascending (deterministic by ref_coord for ties)
    std::vector<const Anchor*> sorted = anchors;
    std::sort(sorted.begin(), sorted.end(),
              [](const Anchor* a, const Anchor* b) {
                  if (a->query_pos == b->query_pos) {
                      // Tie-break by ref_coord for determinism
                      return a->ref_coord < b->ref_coord;
                  }
                  return a->query_pos < b->query_pos;
              });

    std::vector<const Anchor*> top;
    top.reserve(sorted.size());

    std::size_t current_pos = std::numeric_limits<std::size_t>::max();
    for (const auto* a : sorted) {
        if (a->query_pos != current_pos) {
            top.push_back(a);
            current_pos = a->query_pos;
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

// Helper: sample probe anchors with stride
std::vector<const Anchor*> sample_probe_anchors(
    const std::vector<const Anchor*>& anchors,
    std::size_t query_length,
    std::size_t max_probes,
    std::size_t stride_override) {
    if (anchors.empty() || max_probes == 0) return {};

    const std::size_t stride = (stride_override > 0)
                                   ? stride_override
                                   : compute_probe_stride(query_length, max_probes);

    // Reduce to top anchor per query position
    auto top_per_pos = select_top_anchors_per_position(anchors);
    if (top_per_pos.empty()) return {};

    // Sample at stride intervals
    std::vector<const Anchor*> sampled;
    sampled.reserve(std::min(top_per_pos.size(), max_probes));

    std::size_t next_target = 0;
    for (const auto* a : top_per_pos) {
        if (sampled.empty()) {
            sampled.push_back(a);
            next_target = a->query_pos + stride;
            if (sampled.size() >= max_probes) break;
            continue;
        }
        if (a->query_pos >= next_target) {
            sampled.push_back(a);
            next_target = a->query_pos + stride;
            if (sampled.size() >= max_probes) break;
        }
    }

    // Ensure last position is considered if undersampled
    if (!top_per_pos.empty() && sampled.size() < max_probes) {
        const auto* last_anchor = top_per_pos.back();
        if (sampled.empty() || sampled.back()->query_pos != last_anchor->query_pos) {
            sampled.push_back(last_anchor);
        }
    }

    // Cap to max_probes
    if (sampled.size() > max_probes) {
        sampled.resize(max_probes);
    }

    return sampled;
}

// Helper to calculate coverage with interval merging
std::size_t calculate_cluster_coverage(const std::vector<const Anchor*>& cluster) {
    if (cluster.empty()) return 0;

    // Build coverage intervals [start, end]
    std::vector<std::pair<std::size_t, std::size_t>> intervals;
    intervals.reserve(cluster.size());
    for (const auto* a : cluster) {
        std::size_t start = a->query_pos;
        std::size_t end = a->query_pos + a->length - 1;
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
    std::vector<const Anchor*> anchors;  // All anchors in cluster
    std::size_t coverage{0};
    SeedAnchor top_anchor{};
    double top_anchor_score{0.0};
};

// Shared clustering logic (used by both FSE and Probe)
std::vector<ScoredCluster> cluster_and_score_anchors(
    const std::vector<Anchor>& anchors,
    const SeedClustererConfig& config) {

    if (anchors.empty()) return {};

    // 1. Group anchors by path_id
    std::unordered_map<std::size_t, std::vector<const Anchor*>> by_path;
    by_path.reserve(anchors.size());
    for (const auto& a : anchors) {
        by_path[a.path_id].push_back(&a);
    }

    // 2. Within each path, cluster by diagonal with gap tolerance
    std::vector<std::vector<const Anchor*>> all_clusters;

    for (auto& [path_id, path_anchors] : by_path) {
        // Sort by diagonal (GA-style: ref_coord - query_pos)
        std::sort(path_anchors.begin(), path_anchors.end(),
                  [](const Anchor* a, const Anchor* b) {
                      const std::int64_t da = a->ref_coord - static_cast<std::int64_t>(a->query_pos);
                      const std::int64_t db = b->ref_coord - static_cast<std::int64_t>(b->query_pos);
                      if (da == db) return a->query_pos < b->query_pos;  // Tie-break
                      return da < db;
                  });

        // Split clusters when diagonal gap > cutoff
        std::vector<const Anchor*> current;
        current.reserve(path_anchors.size());
        std::int64_t last_diag = 0;
        bool first = true;

        for (const auto* a : path_anchors) {
            const std::int64_t diag = a->ref_coord - static_cast<std::int64_t>(a->query_pos);
            if (first) {
                current.push_back(a);
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
            current.push_back(a);
            last_diag = diag;
        }
        // Add final cluster if meets min size
        if (!current.empty() && current.size() >= config.min_cluster_size) {
            all_clusters.push_back(std::move(current));
        }
    }

    if (all_clusters.empty()) return {};

    // 3. Score anchors in each cluster and find top anchor
    std::vector<ScoredCluster> scored_clusters;
    scored_clusters.reserve(all_clusters.size());

    for (const auto& cluster_anchors : all_clusters) {
        const std::size_t coverage = calculate_cluster_coverage(cluster_anchors);

        // Score based on coverage (anchors don't have frequency info)
        const double score = static_cast<double>(coverage);

        // Select first anchor as representative (they're already sorted by diagonal)
        const auto* best_anchor = cluster_anchors.front();

        SeedAnchor seed_anchor{
            .target = index::SeedHit{
                .node_id = best_anchor->node_id,
                .offset = best_anchor->node_offset,
                .length = best_anchor->length
            },
            .read_pos = best_anchor->query_pos,
            .score = score,
            .cluster_id = 0,  // Will be set later after sorting
        };

        scored_clusters.push_back(ScoredCluster{
            .anchors = cluster_anchors,
            .coverage = coverage,
            .top_anchor = seed_anchor,
            .top_anchor_score = score,
        });
    }

    // 4. Rank clusters by top anchor score
    std::sort(scored_clusters.begin(), scored_clusters.end(),
              [](const ScoredCluster& a, const ScoredCluster& b) {
                  return a.top_anchor_score > b.top_anchor_score;
              });

    // 5. Assign cluster IDs based on rank
    for (std::size_t i = 0; i < scored_clusters.size(); ++i) {
        scored_clusters[i].top_anchor.cluster_id = i;
    }

    // 6. Cap to max_clusters if specified
    if (config.max_clusters > 0 && static_cast<int>(scored_clusters.size()) > config.max_clusters) {
        scored_clusters.resize(static_cast<std::size_t>(config.max_clusters));
    }

    return scored_clusters;
}

class FseClusterer : public AnchorClusterer {
public:
    explicit FseClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<Anchor>& anchors) const override {
        ClusterSummary summary;

        // Anchors already expanded by caller
        summary.expanded_anchor_count = anchors.size();

        // Use shared clustering logic
        auto scored_clusters = cluster_and_score_anchors(anchors, config_);
        if (scored_clusters.empty()) return summary;

        // Return top anchors from selected clusters (multi-mapping support)
        summary.anchors.reserve(scored_clusters.size());
        for (const auto& cluster : scored_clusters) {
            summary.anchors.push_back(cluster.top_anchor);
        }

        // Overall score is the best cluster score
        summary.score = scored_clusters[0].top_anchor_score;

        return summary;
    }

    std::string name() const override { return "fse"; }

private:
    SeedClustererConfig config_;
};

class NoopClusterer : public AnchorClusterer {
public:
    ClusterSummary cluster(const std::vector<Anchor>& anchors) const override {
        ClusterSummary summary;
        summary.expanded_anchor_count = anchors.size();
        summary.anchors.reserve(anchors.size());
        for (const auto& a : anchors) {
            summary.anchors.push_back(SeedAnchor{
                .target = index::SeedHit{
                    .node_id = a.node_id,
                    .offset = a.node_offset,
                    .length = a.length
                },
                .read_pos = a.query_pos,
                .score = 0.0,
            });
        }
        return summary;
    }

    std::string name() const override { return "noop"; }
};

class ProbeClusterer : public AnchorClusterer {
public:
    explicit ProbeClusterer(SeedClustererConfig cfg) : config_(std::move(cfg)) {}

    ClusterSummary cluster(const std::vector<Anchor>& anchors) const override {
        ClusterSummary summary;

        // Anchors already expanded by caller
        summary.expanded_anchor_count = anchors.size();

        // Use shared clustering logic (same as FSE)
        auto scored_clusters = cluster_and_score_anchors(anchors, config_);
        if (scored_clusters.empty()) return summary;

        // Estimate query length from anchors (for probe stride calculation)
        std::size_t query_length = 0;
        for (const auto& a : anchors) {
            query_length = std::max(query_length, a.query_pos + a.length);
        }
        if (query_length == 0) query_length = 1000;  // Fallback if no anchors

        // Sample probes from each cluster and return as grouped clusters
        summary.clusters.reserve(scored_clusters.size());

        for (const auto& scored_cluster : scored_clusters) {
            // Sample probe anchors from this cluster
            auto probes = sample_probe_anchors(
                scored_cluster.anchors,
                query_length,
                config_.max_probes_per_cluster,
                config_.probe_stride);

            if (probes.empty()) continue;

            // Build ClusterGroup with sampled probes
            ClusterGroup group;
            group.cluster_score = scored_cluster.top_anchor_score;
            group.anchors.reserve(probes.size());

            for (const auto* probe : probes) {
                group.anchors.push_back(SeedAnchor{
                    .target = index::SeedHit{
                        .node_id = probe->node_id,
                        .offset = probe->node_offset,
                        .length = probe->length
                    },
                    .read_pos = probe->query_pos,
                    .score = scored_cluster.top_anchor_score,  // Use cluster score
                    .cluster_id = scored_cluster.top_anchor.cluster_id,
                });
            }

            summary.clusters.push_back(std::move(group));
        }

        // Overall score is the best cluster score
        if (!scored_clusters.empty()) {
            summary.score = scored_clusters[0].top_anchor_score;
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

AnchorClustererPtr make_anchor_clusterer(
    const SeedClustererConfig& config,
    const index::GraphStore* graph_store) {

    // Suppress unused parameter warning (kept for potential future use)
    (void)graph_store;

    if (config.backend == "fse" || config.backend.empty()) {
        return std::make_unique<FseClusterer>(config);
    }
    if (config.backend == "probe") {
        return std::make_unique<ProbeClusterer>(config);
    }
    if (config.backend == "dp-chain") {
        // Build DPChainClustererConfig from SeedClustererConfig dp_* fields
        DPChainClustererConfig dp_config;
        dp_config.max_dist = config.dp_max_dist;
        dp_config.max_diag_dev = config.dp_max_diag_dev;
        dp_config.gap_penalty_factor = config.dp_gap_penalty;
        dp_config.diag_penalty_factor = config.dp_diag_penalty;
        dp_config.overlap_penalty_factor = config.dp_overlap_penalty;
        dp_config.anchor_weight = config.dp_anchor_weight;
        dp_config.min_chain_score = config.dp_min_chain_score;
        dp_config.max_chains = config.dp_max_chains;
        LOG_INFO("DP chain config: max_dist=" + std::to_string(dp_config.max_dist) +
                 ", max_diag_dev=" + std::to_string(dp_config.max_diag_dev) +
                 ", gap_penalty=" + std::to_string(dp_config.gap_penalty_factor) +
                 ", diag_penalty=" + std::to_string(dp_config.diag_penalty_factor) +
                 ", overlap_penalty=" + std::to_string(dp_config.overlap_penalty_factor) +
                 ", anchor_weight=" + std::to_string(dp_config.anchor_weight) +
                 ", min_score=" + std::to_string(dp_config.min_chain_score) +
                 ", max_chains=" + std::to_string(dp_config.max_chains));
        return std::make_unique<DPChainClusterer>(dp_config);
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
