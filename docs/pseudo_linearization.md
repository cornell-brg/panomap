# Pseudo-Linearization

## Overview

Pseudo-linearization analyzes the directional `AlnGraph` topology to assign **chain IDs** and **linear coordinates** to nodes. This metadata enables efficient seed clustering (O(n) by chain ID) and colinear chaining (O(n log n) via linear coordinate sorting) during the mapping phase.

The algorithm handles complex graph structures robustly via a multi-step process: SCC detection, tip folding, cycle folding, superbubble detection, and coordinate propagation.

## Goals

1. **Group related nodes**: Assign chain IDs to nodes that belong to the same linear context
2. **Order nodes within chains**: Assign approximate linear positions for colinear chaining
3. **Handle complex topology**: Deal with cycles, tips, and nested variation correctly
4. **Enable efficient mapping**: Provide O(n) clustering and O(n log n) chaining

## Core Concept: Why Linearization?

### The Mapping Problem

During read mapping, we extract seeds (k-mer hits) from the read and graph:

```
Read seeds: [seed1, seed2, seed3, ...]
Graph hits: [(node_5, offset_10), (node_12, offset_3), ...]
```

To form an alignment, we need to:
1. **Cluster seeds**: Group seeds that likely come from the same alignment
2. **Chain seeds**: Find colinear sets of seeds with consistent spacing

### Without Linearization: O(n²) Clustering

```
For each pair of seeds (i, j):
    if graph_distance(seed_i.node, seed_j.node) < threshold:
        cluster together
```

This requires computing all-pairs graph distances: **O(n²)** per read!

### With Linearization: O(n) Clustering + O(n log n) Chaining

```
// O(n) clustering by chain ID
for seed in seeds:
    clusters[seed.chain_id].append(seed)

// O(n log n) chaining per cluster (sort by linear coordinate)
for cluster in clusters:
    cluster.sort(key=lambda s: s.linear_coord)
    chain_with_gap_costs(cluster)
```

Linear coordinates provide a **total ordering** within each chain, enabling efficient sorting and gap-cost computation without graph traversal.

## Algorithm Pipeline

Six stages run in order, using union-find to group nodes into chains:

```
AlnGraph
  → SCC Detection (computeScc)
  → Tip Folding (chainTips)
  → Cycle Folding (chainCycles)
  → Superbubble Detection (chainSuperbubbles)
  → Chain ID Assignment
  → Linear Coordinate Assignment (assignLinearPositions)
```

---

## Step 1: SCC Detection (`computeScc`)

### Purpose

Identify **strongly connected components** (SCCs) — regions of the graph where all nodes are mutually reachable via directed paths. SCCs represent cyclic structures (tandem repeats, complex SVs) that cannot be linearized as simple DAG chains.

### Why This is Needed

The superbubble algorithm (Step 4) assumes **DAG structure**. Starting detection from nodes in cycles will:
- Produce incorrect results (false superbubbles)
- Cause infinite loops (cycles have no topological ordering)

By identifying SCCs upfront, we can **skip cyclic nodes** during superbubble detection.

Additionally, SCC detection provides **component numbers** — a topological ordering where downstream components have higher numbers. This is used in tip detection and superbubble chaining to ensure we only chain nodes in compatible topological layers.

### Algorithm: Iterative Tarjan's SCC

**Tarjan's algorithm** finds SCCs in O(V+E) time via DFS with low-link values.

**Key idea**:
- Track `index[v]`: Discovery time of node v
- Track `lowlink[v]`: Lowest index reachable from v via back-edges
- If `lowlink[v] == index[v]`, then v is the **root** of an SCC

**Implementation details**:

1. **Iterative (not recursive)**: explicit call stack to avoid recursion limits
2. **State machine**: track node + next child to visit
3. **Stack-based SCC collection**: when root found, pop stack to build component

**Example**:

```
Graph: A → B → C → A (cycle)
       D → A

Tarjan discovers:
  Component 0: {A, B, C}  (nontrivial SCC, size > 1)
  Component 1: {D}        (trivial SCC, size = 1)

After reversal:
  Component 0: {D}        (upstream)
  Component 1: {A, B, C}  (downstream)
```

### Component Number Reversal

**Critical step**: Tarjan assigns components in **reverse topological order** (later discoveries get lower numbers). We reverse this to get **forward topological order**:

```cpp
for (v : nodes) {
    component[v] = (total_components - 1) - component[v];
}
```

**Why**: Tip folding and superbubble chaining rely on the invariant that edges go from **lower** to **higher** component numbers:

```cpp
if (component[end] <= component[start]) {
    skip;  // Not a valid forward bubble
}
```

### Nontrivial SCC Marking

An SCC is **nontrivial** if:
- Size > 1, OR
- Size = 1 AND has a self-loop

Nodes in nontrivial SCCs are marked `in_scc[v] = true` and **skipped** in superbubble detection.

### Output

```cpp
struct SccResult {
    std::vector<std::size_t> component;  // component ID per node (topologically ordered)
    std::vector<bool> in_scc;            // true if in nontrivial SCC or self-loop
    std::size_t components;              // total component count
};
```

### Validation

After reversal, verify topological ordering:

```cpp
for edge (v → u):
    assert(component[u] >= component[v]);  // Edges go forward
```

If this fails, the SCC algorithm has a bug.

---

## Step 2: Tip Folding (`chainTips`)

### Purpose

Identify and chain **tip components** — dead-end branches that have no continuation beyond their component boundary.

### What is a Tip?

A **forward tip** is a component with:
- No outgoing edges to non-tip components in the same or downstream components

A **backward tip** is a component with:
- No incoming edges from non-tip components in the same or upstream components

**Example**:

```
Main path: A → B → C
Tip:       D → E (ends, no continuation)

E is a forward tip (no outgoing edges beyond component boundary)
```

### Why Fold Tips?

Tips represent:
- Assembly artifacts (dead-end contigs)
- Singleton variants (private mutations with no downstream context)
- Incomplete graph regions

They don't participate in variation bubbles (which require reconvergence). Folding them early prevents interference with superbubble detection.

### Algorithm

1. **Forward tip detection** (reverse component order): a component is a forward tip if all outgoing edges stay within components already marked as tips (including itself).
2. **Backward tip detection** (forward component order): a component is a backward tip if all incoming edges originate from components already marked as tips (including itself).
3. **Union-find merging**: union edges within tip components; mark `ignorable_tip[v] = true` for tip nodes.

### Output

- **Union-find structure** (`parent`, `rank`) with tip nodes merged
- **`ignorable_tip[v]`**: Boolean flag for nodes in tips (used to skip them in superbubble detection)

### Why Component Order Matters

Traversing in **reverse order** for forward tips ensures we detect tips **leaf-first**:

```
Components: 0 → 1 → 2 → 3 (topological order)
            A   B   C   D

Traverse reverse: 3, 2, 1, 0

Step 1: Check D (component 3)
  - No outgoing edges → forward tip
  - Mark component 3 as tip

Step 2: Check C (component 2)
  - Outgoing edges only to D (already marked tip) → forward tip
  - Mark component 2 as tip

Step 3: Check B (component 1)
  - Outgoing edges to C (tip) and A (non-tip) → NOT a tip
```

---

## Step 3: Cycle Folding (`chainCycles`)

### Purpose

Detect and collapse **simple 2-cycles** (bidirectional edges between two nodes).

### What is a 2-Cycle?

Two nodes `v` and `u` form a 2-cycle if:
- `v` has **unique** forward neighbor `u` (ignoring tips and self-loops)
- `u` has **unique** backward neighbor `v`

**Example**:

```
v ⇄ u  (bidirectional edges)
```

### Why Fold Cycles?

Simple cycles represent:
- Bidirectional structural variants (inversions)
- Assembly artifacts (conflicting orientations)

Collapsing them simplifies the graph for superbubble detection by removing cyclic dependencies.

### Algorithm

For each node `v`:
1. Find unique forward neighbor `u_fw` (skip tips, self-loops, multi-neighbors)
2. Find unique backward neighbor `u_bw`
3. If `u_fw == u_bw` and both are unique: union `v` with `u_fw` and mark `v` ignorable.

### Output

- Updated union-find structure
- Updated `ignorable_tip` flags

---

## Step 4: Superbubble Detection and Chaining (`chainSuperbubbles`)

### Purpose

Identify **superbubbles** — maximal DAG regions with a single entry and single exit node — and chain all nodes within each bubble.

### What is a Superbubble?

A superbubble is a pair of nodes `(s, t)` where:
- `s` has out-degree > 0
- All paths from `s` reconverge at `t` before reaching any other node
- No cycles involve `s` or nodes between `s` and `t`

**Example**:

```
     ┌→ B ─┐
  s ─┤     ├→ t
     └→ C ─┘

(s, t) is a superbubble
```

Superbubbles represent **variation regions** (SNPs, indels, SVs) where multiple alleles diverge and reconverge.

### Algorithm

For each node `s` (excluding SCC nodes and ignorable tips):

1. Skip if `in_scc[s]` is true.
2. Run superbubble detection from `s` (Onodera et al. 2013), skipping ignorable tips, return exit `t` if found.
3. Require `component[t] > component[s]` (strictly downstream).
4. Union `s`/`t`, traverse interior (reachable from `s` without passing `t`), and union all interior with `s`.

### Superbubble Detection Algorithm (Onodera et al. 2013)

**Key idea**: A superbubble `(s, t)` satisfies:
- `t` is the first node where all paths from `s` reconverge

**Implementation** (simplified):

```
visited = {s}
reachable = {successors of s}

while reachable not empty:
    candidate_t = node in reachable with no unvisited predecessors

    if candidate_t has predecessors outside visited:
        return None  // Not a valid superbubble

    if all paths reconverge at candidate_t:
        return candidate_t  // Found exit node

    visited.add(candidate_t)
    reachable.update(successors of candidate_t)
```

### Why Skip Tips?

Ignorable tips are excluded from superbubble detection to avoid false positives:

```
     ┌→ B ─┐
  s ─┤     ├→ t
     └→ D (tip, dead-end)

Without skipping tip: Algorithm might incorrectly identify (s, B) as a bubble
With skipping tip: Correctly identifies (s, t) as the superbubble
```

### Output

- Updated union-find structure with superbubble nodes merged into chains

---

## Step 5: Chain ID Assignment (`assignChainIds`)

### Purpose

Convert the union-find structure into explicit, **deterministic** chain IDs. This provides a dense, sequential numbering (0, 1, 2, ...) for downstream processing and ensures reproducible results across runs.

### Algorithm

**Implementation:** `assignChainIds(const UnionFind& uf)` (pseudo_linearize.cpp:475-502)

```cpp
1. Path-compress union-find to get representative for each node:
   for (i in 0..n):
       reps[i] = uf.find(i)

2. Build sorted unique list of representatives:
   unique_reps = sorted(unique(reps))

3. Map each representative to a sequential chain ID:
   for (cid in 0..unique_reps.size()):
       rep_to_chain[unique_reps[cid]] = cid

4. Assign each node the chain ID of its representative:
   for (i in 0..n):
       chain_id[i] = rep_to_chain[reps[i]]
```

### Why Sort Representatives?

**Determinism**: Without sorting, chain IDs would depend on the order representatives are encountered (which depends on node iteration order and union-find structure). Sorting by representative node ID ensures:

- Same graph → same chain IDs
- Reproducible across runs
- Easier debugging (lower node IDs → lower chain IDs)

**Legacy compatibility**: The legacy PIRU implementation (`chain_superbubbles_ga`, lines 717-724) sorts representatives, and we match this behavior exactly.

### Example

```
Union-find after superbubble chaining:
  Nodes 0, 3, 7 → representative 0
  Nodes 1, 5 → representative 1
  Nodes 2, 4, 6 → representative 2

Representatives: [0, 1, 2]  (already sorted)

Chain ID assignment:
  rep 0 → chain 0
  rep 1 → chain 1
  rep 2 → chain 2

Result:
  chain_id[0] = 0, chain_id[3] = 0, chain_id[7] = 0
  chain_id[1] = 1, chain_id[5] = 1
  chain_id[2] = 2, chain_id[4] = 2, chain_id[6] = 2
```

### Output

- **Return value**: `std::vector<std::size_t> chain_id` (size = num_nodes)
- **Range**: 0 to num_chains-1 (dense, sequential)
- **All nodes receive a chain ID** (singleton chains for nodes not merged)

---

## Step 6: Linear Coordinate Assignment (`assignLinearPositions`)

### Purpose

Assign approximate linear positions within each chain to enable colinear chaining during mapping.

### Algorithm

**Implementation:** `assignLinearPositions(graph, chain_ids, scc)` (pseudo_linearize.cpp:504-575)

**Signature:**
```cpp
std::vector<std::int64_t> assignLinearPositions(
    const AlnGraph& graph,
    const std::vector<std::size_t>& chain_ids,
    const SccResult& scc);
```

**Steps:**

```cpp
1. Build chain membership:
   for (v in 0..n):
       chain_members[chain_ids[v]].push_back(v)

2. Initialize all positions to -1 (unassigned)

3. For each chain:
   a. Choose root node for traversal (component-aware):
      - If SCC data available: pick node with lowest (component_number, node_id)
      - Otherwise: pick node with lowest node_id

   b. Initialize DFS stack:
      - positions[root] = INT32_MAX (2^31 - 1)
      - stack.push((root, INT32_MAX))

   c. DFS traversal (stack-based):
      while stack not empty:
          (v, dist) = stack.pop()

          if positions[v] != -1: continue  // Already visited
          positions[v] = dist

          // Forward edges: v → u
          for u in graph.outgoing(v):
              if chain_ids[u] != chain_ids[v]: continue
              if positions[u] != -1: continue
              stack.push((u, dist + len(v)))

          // Backward edges: u → v
          for u in graph.incoming(v):
              if chain_ids[u] != chain_ids[v]: continue
              if positions[u] != -1: continue
              stack.push((u, dist - len(u)))
```

### Component-Aware Root Selection

**Why?** Choosing the root node based on component number (from SCC detection) provides better starting points for propagation:

- **Topologically earlier nodes** are preferred (start from "upstream" in the graph)
- **Consistent with graph flow** (edges go from low→high component numbers)
- **Deterministic** tie-breaking via node ID

**Algorithm (lines 527-541):**
```cpp
root = members[0]
for v in members:
    if scc.components > 0:
        if (component[v] < component[root]) OR
           (component[v] == component[root] AND v < root):
            root = v
    else:
        if v < root:
            root = v
```

**Example:**
```
Chain members: {5, 12, 3, 8}
Component numbers: {1, 0, 1, 0}

Comparison:
  v=12: component=0, root=5: component=1 → 0 < 1 → root=12
  v=3:  component=1, root=12: component=0 → 1 > 0 → keep root=12
  v=8:  component=0, root=12: component=0 → tie, 8 < 12 → root=8

Final root: 8 (lowest component number, then lowest node ID)
```

### Stack-Based DFS Traversal

**Why DFS (not BFS)?**
- **Matches legacy implementation** (uses vector with back/pop_back)
- **Memory efficient** for deep chains
- **Order doesn't affect correctness** (all nodes reached)

**Visited tracking:** Uses `positions[v] != -1` instead of separate visited set (saves memory)

### Properties

- **Monotonic along paths**: Nodes along the same path have increasing positions
- **Overlapping in bubbles**: Nodes in parallel bubble arms may have overlapping positions
- **Distance approximation**: Position difference approximates sequence length along edges
- **May be negative**: Backward propagation can assign negative positions (hence int64_t)

### Example

```
Chain: A → B → C
       A → D → C  (bubble)

Seed: A (position = 100000)

Forward propagation:
  B: 100000 + len(A) = 100050
  D: 100000 + len(A) = 100050  (same as B, parallel path)
  C: 100050 + len(B) = 100100  (from B)
  C: 100050 + len(D) = 100120  (from D, overwrites previous)

Result:
  A: 100000
  B: 100050
  D: 100050  (overlaps with B)
  C: 100120
```

### Why Approximate Positions?

Linear coordinates don't need to be **exact** — they just need to provide a **total ordering** for sorting:

```
Read seeds: [(node_B, offset_5), (node_D, offset_10), (node_C, offset_3)]

Sort by linear_coord + offset:
  (B, 5): 100050 + 5 = 100055
  (D, 10): 100050 + 10 = 100060
  (C, 3): 100120 + 3 = 100123

Sorted order: B, D, C → enables gap-cost computation
```

### Output

- **Return value**: `std::vector<std::int64_t> positions` (size = num_nodes)
- **Value range**: Can be positive, zero, or negative (backward propagation)
- **Unassigned nodes**: Position = -1 (if chain has issues, should not occur normally)

---

## Final Output

After all 6 steps, the pseudo-linearization pipeline produces two vectors:

```cpp
// From assignChainIds()
std::vector<std::size_t> chain_ids;  // Size = num_nodes, range [0, num_chains)

// From assignLinearPositions()
std::vector<std::int64_t> positions;  // Size = num_nodes, arbitrary int64_t values
```

These vectors are **parallel arrays** indexed by node ID:
```cpp
std::size_t node_id = 42;
std::size_t chain = chain_ids[node_id];      // Which chain does this node belong to?
std::int64_t pos = positions[node_id];       // What's its approximate linear position?
```

### Storage in AlnGraph

The `AlnNode` struct stores these values as optional fields:

```cpp
struct AlnNode {
    std::size_t id;
    std::string label;
    std::string original_id;
    bool is_reverse;
    std::string sequence;
    std::optional<std::int64_t> chain_id;       // Set after pseudo-linearization
    std::optional<std::int64_t> linear_position; // Set after pseudo-linearization
};
```

**Assignment:**
```cpp
SccResult scc = computeScc(graph);
TipFoldingResult tips = chainTips(graph, scc);
chainCycles(graph, tips);
SuperbubbleResult sb = chainSuperbubbles(graph, scc, tips);
auto chain_ids = assignChainIds(sb.uf);
auto positions = assignLinearPositions(graph, chain_ids, scc);

// Store in graph
for (std::size_t i = 0; i < graph.nodeCount(); ++i) {
    graph.node(i).chain_id = chain_ids[i];
    graph.node(i).linear_position = positions[i];
}
```

**Properties:**
- **All nodes receive a chain ID** (singleton chains for un-merged nodes)
- **All nodes receive a position** (except on error: position = -1)

---

## Usage During Mapping

### Seed Clustering (O(n))

```cpp
std::unordered_map<ChainId, std::vector<Seed>> clusters;
for (seed : seeds) {
    ChainId cid = graph_store.chainId(seed.node_id);
    clusters[cid].push_back(seed);
}
```

### Colinear Chaining (O(n log n) per cluster)

```cpp
for (auto& [cid, seeds] : clusters) {
    // Sort seeds by linear position
    std::sort(seeds.begin(), seeds.end(), [&](auto& a, auto& b) {
        return graph_store.linearPosition(a.node_id) <
               graph_store.linearPosition(b.node_id);
    });

    // Apply gap-cost chaining
    chain_with_gaps(seeds);
}
```

---

## Biological Motivation

### SCCs: Tandem Repeats

```
Graph: A → B → C → B → C → B → C → D
            └────┘   └────┘   └────┘

SCC: {B, C} (nontrivial, size = 2)
```

Tandem repeats cycle indefinitely — no single linear traversal order exists. The SCC detection correctly identifies them as non-linearizable.

### Tips: Private Variants

```
Main path: A → B → C
Tip:       X → Y (only in one sample, dead-end)
```

Private variants that don't rejoin the main graph are folded as tips.

### Superbubbles: Variation Regions

```
     ┌→ REF allele ─┐
  s ─┤              ├→ t
     └→ ALT allele ─┘
```

SNPs, indels, and SVs create superbubbles. Chaining them enables colinear alignment within variation contexts.

---

## Implementation Details

### Union-Find with Path Compression

```cpp
std::size_t find(std::size_t x, std::vector<std::size_t>& parent) {
    if (parent[x] != x) {
        parent[x] = find(parent[x], parent);  // Path compression
    }
    return parent[x];
}

void union_sets(std::size_t x, std::size_t y,
                std::vector<std::size_t>& parent,
                std::vector<std::size_t>& rank) {
    x = find(x, parent);
    y = find(y, parent);
    if (x == y) return;

    // Union by rank
    if (rank[x] < rank[y]) {
        parent[x] = y;
    } else if (rank[x] > rank[y]) {
        parent[y] = x;
    } else {
        parent[y] = x;
        rank[x]++;
    }
}
```

### Topological Ordering Invariant

Throughout the pipeline, component numbers must satisfy:

```
For all edges (v → u): component[v] <= component[u]
```

This invariant is validated after SCC detection and relied upon in tip folding and superbubble chaining.

---

## Performance

- **Time complexity**: O(V + E) for each step → O(V + E) overall
- **Space complexity**: O(V) for auxiliary arrays (index, lowlink, stack, union-find)

For a graph with 10M nodes and 20M edges, pseudo-linearization completes in ~10 seconds.

---

## References

- **Implementation**: `piru/include/index/pseudo_linearize.hpp`, `piru/src/index/pseudo_linearize.cpp`
- **Legacy reference**: `PIRU-workspace/PIRU/src/graph/alignment_graph.cpp` (lines 178-871)
- **Design doc**: `piru/docs/graph_indexing.md` (lines 139-235)
- **Tarjan's algorithm**: Tarjan, Robert (1972). "Depth-first search and linear graph algorithms"
- **Superbubble algorithm**: Onodera et al. (2013). "Defining and identifying superbubbles" ([arXiv:1307.7925](https://arxiv.org/pdf/1307.7925))
