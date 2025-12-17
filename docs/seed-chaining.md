# Seed Chaining Pipeline

Path-guided seed chaining turns raw SeedStore hits into colinear anchor chains ready for alignment extension. The implementation linearizes nodes along reference paths, expands seeds to anchors, optionally merges overlapping anchors, and uses DP chaining to pick the best chain(s) per read.

## Pipeline at a Glance

```
  SeedStore          PathIndex
    Lookup         (linearization coords)
      │                 │
      ▼                 │
  Seed Hits             │
(q, v, off, l)          │
      │                 │
      └────────┬────────┘
               ▼
      ┌───────────────┐
      │    Anchor     │
      │  Expansion    │
      └───────────────┘
               │
               ▼
          Anchors
      (q, r, l, p)
      [one seed → many anchors]
               │
               ▼
      ┌───────────────┐
      │    Anchor     │
      │   Merging     │
      │ (same path)   │
      └───────────────┘
               │
               ▼
       Merged Anchors
      (q, r, l', p)
               │
               ▼
      ┌───────────────┐
      │  DP Chaining  │
      │  (colinear)   │
      └───────────────┘
               │
               ▼
       Anchor Chains
     (multiple chains
      per read for
      multi-mapping)
```

## Core Data

- **SeedHit** `(query_pos, node_id, node_offset, length)` from SeedStore lookup.
- **Path occurrence index** `(node_id, offset) → [(path_id, ref_coord)]` built during indexing by walking reference paths and recording cumulative coordinates.
- **Anchor** `(query_pos, ref_coord, length, path_id [, node_id, node_offset])` produced by expansion; one seed can emit many anchors (one per path occurrence).

## Stage 1: Anchor Expansion

Use the occurrence index to map each seed to linear coordinates:

```
for seed in seed_hits:
  for (path_id, ref_coord) in occurrence_index[(seed.node_id, seed.node_offset)]:
    anchors.push({
      query_pos: seed.query_pos,
      ref_coord: ref_coord,
      length: seed.length,
      path_id: path_id,
      node_id: seed.node_id,
      node_offset: seed.node_offset
    })
```

Nodes on multiple haplotypes or repeated on a path produce multiple anchors, enabling haplotype-aware chaining and multi-mapping.

## Stage 2: (Optional) Anchor Merging

Goal: collapse adjacent or overlapping anchors on the same path to shrink the workload and reflect contiguous similarity.

Process:
- Group by `path_id`, sort by `(ref_coord, query_pos)`.
- Merge when both query and ref gaps are within tolerance (default 0).
- Reject merges if either axis jumps backwards beyond tolerance.

Result: fewer, longer anchors with updated `length`.

## Stage 3: Sorting

Sort anchors by `(path_id, ref_coord, query_pos)` so DP traversal sees increasing reference order per path (and deterministic tie-breaking).

## Stage 4: DP Colinear Chaining

Objective: pick colinear anchors that maximize chain score.

Constraints per candidate predecessor `j` → `i`:
- Order: `q_j < q_i` and `r_j < r_i`
- Distance: `dq = q_i - q_j ≤ max_dist`, `dr = r_i - r_j ≤ max_dist`
- Diagonal: `|dq - dr| ≤ max_diag_dev`
- Path: same path unless cross-haplotype chaining is enabled (path switch penalty applied)

Scoring (per edge `j → i`):
```
score = dp[j] + anchor_weight * l_i
        - gap_cost(dq, dr, l_j, l_i)
        - (p_j != p_i ? path_switch_cost : 0)
```

Gap cost penalizes distance, off-diagonal drift, and overlap:
```
query_gap = dq - l_j
ref_gap   = dr - l_j
avg_gap   = avg(max(0, query_gap), max(0, ref_gap))
diag_dev  = abs(dr - dq)
avg_ovlp  = avg(max(0, -query_gap), max(0, -ref_gap))

gap_cost = gap_penalty_factor * avg_gap
         + diag_penalty_factor * diag_dev
         + overlap_penalty_factor * avg_ovlp
```

Chain extraction:
- Run DP once over sorted anchors.
- Repeatedly pick the best unused endpoint (skipping intervals already covered on the same path with a small margin), backtrack via `pred`, and mark anchors used until `max_chains` or `min_chain_score` stops the loop.
- Chains may share prefixes; same-path redundancy is filtered by interval coverage.

## Key Parameters (defaults)

- `merge_tolerance = 0`
- `max_dist = 5000`
- `max_diag_dev = 500`
- `allow_cross_haplotypes = false`
- `path_switch_cost = 50.0`
- `min_chain_score = 100`
- `max_chains = 10`
- `anchor_weight = 1.0`
- `gap_penalty_factor = 0.1`
- `diag_penalty_factor = 0.5`
- `overlap_penalty_factor = 2.0`

## Implementation Notes

- Expansion: `src/mapping/anchor_expander.cpp`
- Optional merging: `src/mapping/anchor_merger.cpp`, `include/mapping/anchor_merger.hpp`
- DP chaining: `src/mapping/dp_chain_clusterer.cpp`
- Wiring in mapping pipeline (lookup → expand → merge → chain): `src/mapping/batch_mapper.cpp`

Useful doc links: `docs/pseudo_linearization.md` (superbubble linearization) and `docs/graph_indexing.md` (path indexing).
