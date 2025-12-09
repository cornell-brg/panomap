# Graph Indexing Pipeline

This document details the architecture for the `piru index` command, which processes a variation or assembly graph into a set of indexed data structures (`GraphStore`, `SignalStore`, `SeedStore`) used for mapping.

## High-Level Flow

The diagram below illustrates the major stages of the indexing pipeline, from input graph to final on-disk index components.

**Note**: The pipeline supports two graph types with different transformation strategies:
- **Variation Graph (VG)**: No inherent overlaps between nodes → requires VG→DBG expansion to generate k-1 overlaps
- **De Bruijn Graph (DBG)**: Inherent k-1 overlaps already present → direct transformation preserving overlaps

```
   GFA/vg          pore kmer
   graph            model
     │               │
     └───────┬───────┘
             ▼
      ┌─────────────┐
      │ GraphLoader │
      └─────────────┘
             │
             ▼
       ImportedGraph
      (bidirectional)
             │
             │
    ┌────────┴────────┐
    │                 │
    ▼                 ▼
┌─────────┐       ┌─────────┐
│ VG Path │       │DBG Path │
└─────────┘       └─────────┘
    │                 │
    ▼                 ▼
┌──────────────────────────┐ ┌──────────────────────────┐
│ VG→DBG Transformation    │ │ DBG Transformation       │
│  • Node splitting (F/R)  │ │  • Node splitting (F/R)  │
│  • Internode traversal   │ │  • Preserve k-1 overlaps │
│  • k-1 overlap generation│ │  • No expansion needed   │
└──────────────────────────┘ └──────────────────────────┘
    │                 │
    └────────┬────────┘
             ▼
         AlnGraph
 (directional, overlapping)
             │
             ▼
 ┌──────────────────────────┐
 │ Pseudo-Linearization     │
 │  • Superbubble detection │
 │  • Chain ID assignment   │
 │  • Linear coord mapping  │
 └──────────────────────────┘
             │
             ▼
         AlnGraph
  (with chain IDs and
   linear coordinates)
             │
             ▼
 ┌──────────────────────────┐
 │ Squigglize               │
 │ (base sequence → raw     │
 │  signal samples)         │
 └──────────────────────────┘
             │
             ▼
      ┌──────────────┐
      │ Quantization │
      └──────────────┘
             │
             ├─────────────────┐
             │                 │
             ▼                 ▼
     Fuzzy Quantized     Alignment Quantized
       Samples               Samples
     (for seeding)         (for storage)
             │                 │
             ▼                 │
    ┌────────────────┐         │
    │ Seed Generation│         │
    │  & Indexing    │         │
    └────────────────┘         │
             │                 │
             │ produces three components:
             │
   ┌─────────┼─────────────────┤
   │         │                 │
   ▼         ▼                 ▼
SeedStore GraphStore         SignalStore
(hash:     (topology)         (per-node
 seed→                      aln-quant
 node+pos)                   samples)
   │         │                 │
   └─────────┴─────────────────┘
             │
             ▼
          index
       (on disk)
```

## Pipeline Stages & Components

### Input Layer

- **Inputs**:
    1.  **Graph File**: A graph in `GFA` or `vg` format.
    2.  **Graph Type**: A parameter (`--graph-type`) specifying the graph's nature, either `vg` (variation graph) or `dbg` (de Bruijn graph). This dictates the transformation strategy.
    3.  **Pore Model**: An optional `.model` file. If not provided, a built-in model (e.g., R9.4) is used. The model's k-mer size (`k`) is critical for the transformation process.

### Graph Transformation (`ImportedGraph` -> `AlnGraph`)

The initial `ImportedGraph` loaded from disk is bidirectional and may not have the sequence overlap required for signal processing. This stage transforms it into a directional `AlnGraph` with properties suitable for squigglization and mapping.

**Squigglization Context Requirement**: To squigglize (convert base sequences to signal), we slide a k-mer window of size `k` (from the pore model) across each node's sequence and consult the pore model lookup table. At node boundaries, the sliding window needs access to the next `k-1` bases from successor nodes. Without overlaps, we would need to traverse edges and fetch internode sequences during squigglization—an expensive and complex operation. Instead, we precompute overlaps during transformation.

The pipeline supports **two graph types** with different transformation strategies:

#### 1. De Bruijn Graph (DBG) Transformation

DBGs inherently encode overlaps between nodes: by construction, adjacent nodes in a DBG share `k-1` bases. This property makes them naturally suited for squigglization.

-   **Node Splitting (F/R)**: Each bidirectional node is split into two directional nodes (forward and reverse complement). Edges are duplicated and reconnected to maintain graph topology.
-   **Preserve k-1 Overlaps**: The transformation preserves the inherent `k-1` overlaps already present in the DBG structure. Node sequences retain their original overlap content.
-   **No Expansion Needed**: Since overlaps are already present, no additional sequence traversal or node expansion is required.

#### 2. Variation Graph (VG) → DBG Transformation

Variation graphs typically do not have overlapping nodes—adjacent nodes represent distinct sequence regions without shared bases. This poses a problem for squigglization, which requires `k-1` bases of context at node boundaries.

**Transformation Strategy: Path-Guided with Local Expansion Fallback**

VG graphs often contain embedded paths (haplotypes) representing known biological sequences. The implemented transformation leverages these paths for accurate context while providing a fallback mechanism for nodes not covered by any path.

**Implementation**: `PathGuidedTransform` in `piru/include/index/path_guided_transform.hpp`

**Algorithm Overview:**

1. **Bidirectional → Directional Split**: Each bidirectional node is split into forward and reverse complement variants (same as DBG transformation).

2. **Path Walking with Context Addition**: For each embedded path in the VG graph:
   - Trace the path through the directional graph
   - For each step, append `k-1` bases from the next node in the path
   - Create AlnNode with sequence = `base_sequence + next_node[0:k-1]`

3. **Node Duplication by Path Context**: When the same node appears in multiple paths:
   - If successor context is identical → reuse the same AlnNode (share across paths)
   - If successor context differs → create a new AlnNode variant with the different context
   - Track variants via key: `(original_id, successor_context, is_reverse)`

4. **Coverage Tracking**: Track which original nodes were visited by any path. Nodes not covered require fallback handling.

5. **Local Expansion Fallback**: For uncovered nodes (not visited by any path):
   - Perform greedy depth-limited traversal to collect `k-1` context from successors
   - Create one variant per immediate successor
   - Follow first available path to avoid exponential explosion in highly branching regions

**Example**:

```
Original VG with paths:
    Node A (GATTACA)
      ├─→ B (CCC)
      └─→ C (GGG)

    Path P1: A → B
    Path P2: A → C

After Path-Guided Transformation (k=5, k-1=4):
    A_B (GATTACA + CCCC) [from path P1] ──→ B (CCC)
    A_C (GATTACA + GGGG) [from path P2] ──→ C (GGG)

    Two variants of A, each with biologically accurate context from real haplotypes
```

**Example with Node Sharing**:

```
Original VG:
    Node X (AAAA) → Node Y (TTTT) → Node Z (GGGG)

    Path P1: X → Y → Z
    Path P2: X → Y → Z  (same successor contexts)

After Transformation:
    X_Y (AAAA + TTTT) [shared by P1, P2] ──→ Y_Z (TTTT + GGGG)

    Only one variant of each node created, shared across both paths
```

**Handling Uncovered Nodes**:

Nodes not visited by any embedded path use local expansion:

```
Uncovered node N with successors S1, S2:
  → Creates N_S1 with context from S1 (following first path to depth k-1)
  → Creates N_S2 with context from S2 (following first path to depth k-1)
```

**Graph Size Impact**:

- **Best case**: All nodes covered by paths with identical contexts → minimal duplication
- **Typical case**: Node duplication proportional to number of distinct contexts in embedded paths
  - Example: drb1.vg with 12 HLA haplotype paths → node may appear up to 12 times if all paths provide different successor contexts
  - Actual: drb1.vg: 5111 original nodes → 15333 transformed nodes (3x expansion)
- **Worst case**: Highly branching uncovered regions → expansion based on successor count

**Trade-offs**:

- ✅ **Biologically accurate**: Context comes from real haplotype sequences
- ✅ **Efficient for path-covered nodes**: Shares nodes when paths agree on context
- ✅ **Complete coverage**: Local expansion ensures all nodes get valid context
- ❌ **Graph size increase**: More nodes than original VG (but less than naive splitting)
- ✅ **Simple squigglization**: No runtime edge traversal needed

**Handling Special Cases**:

- **Nodes with no successors (tips)**: Get variant with empty context (base sequence only)
- **Nodes in cycles**: Visited set tracking prevents infinite loops during expansion
- **High branching factor**: Greedy traversal (first available path) prevents exponential explosion
- **N bases in sequences**: Emit NaN sentinels during squigglization; seed extraction skips these windows

**Output**: A DBG-like `AlnGraph` where all nodes are directional and have sequences that include `k-1` bases of overlap from their successors, suitable for squigglization. Nodes are duplicated only when their path contexts genuinely differ, preserving biological accuracy.

#### Handling N Bases (Unknown Sequences)

VG graphs often contain `N` bases representing unknown or ambiguous sequence regions. Since the pore model only defines signals for ACGT k-mers, k-mers containing N bases cannot be squigglized into valid signals. The pipeline uses a **sentinel value approach** to handle this:

**Sentinel Value Strategy:**

1. **Squigglization** (`squigglize.cpp`):
   - Check each k-mer for N bases before model lookup
   - If N detected: emit `std::numeric_limits<float>::quiet_NaN()` instead of a signal value
   - Maintains 1:1 position correspondence: `signal[i]` always maps to sequence position `i`

2. **Quantization**:
   - **Float stages** (normalized signal, float alignment): Pass NaN through unchanged
   - **Int stages** (fuzzy int16, alignment int16/int8): Map NaN → `std::numeric_limits<T>::min()`
     - `int16_t`: NaN → `-32768` (min value reserved as sentinel)
     - `int8_t`: NaN → `-128` (min value reserved as sentinel)
     - Valid values clamped to `[min()+1, max()]` to avoid collision

3. **Seed Extraction** (`kmer_seed_extractor.cpp`):
   - Before generating seed hash, check if window contains any sentinel values
   - Skip seed extraction for windows overlapping sentinel positions
   - Result: **Only valid ACGT k-mers enter the seed index**

**Example**: drb1.vg (HLA graph with many N bases)
- 5111 original nodes with sequences containing N bases
- Squigglization emits NaN at N positions (no warnings)
- Seed index contains only 45 unique valid seeds (N-containing windows skipped)
- Clean indexing with INFO-level logging: "Node X: Y k-mers with N bases (marked as NaN)"

**Semantic Meaning**: Unknown base → missing signal (gap) → no seed → these regions are effectively excluded from seeding but preserved in topology.

**Alternative Future Approach**: Tree-based model approximation (deferred):
- For ambiguous k-mers like `AACTN`, compute average of all variants: `avg(model[AACTG], model[AACTA], model[AACTC], model[AACTT])`
- Mark with quality scores indicating uncertainty
- Allows approximate matching in N-rich regions (more complex, requires model changes)

### Pseudo-Linearization

This stage analyzes the directional `AlnGraph` topology to identify linear structures (chains) and assign metadata that enables efficient seed clustering and colinear chaining during mapping. The algorithm follows a multi-step process to handle complex graph structures robustly.

#### Algorithm Overview

The pseudo-linearization follows this pipeline:

1. **SCC (Strongly Connected Components) Detection** via iterative Tarjan's algorithm
2. **Tip Folding** (forward and backward tip detection and chaining)
3. **Simple Cycle Folding** (2-cycle detection and collapsing)
4. **Superbubble Detection and Chaining** (DAG regions only)
5. **Linear Coordinate Assignment** via breadth-first propagation

#### Step 1: SCC Detection (`computeScc`)

-   **Purpose**: Pre-filter cyclic regions and establish topological ordering for subsequent chaining steps.
-   **Algorithm**: Iterative Tarjan's algorithm (stack-based to avoid recursion depth issues on large graphs).
-   **Output**:
    -   `node_component_number[v]`: Component ID for each node, topologically ordered (higher component numbers for downstream nodes)
    -   `node_in_scc[v]`: Boolean flag indicating if node is in a nontrivial SCC (size > 1) or has a self-loop
-   **Why This is Needed**:
    -   **Cycles break linearization**: In cyclic regions (e.g., `A → B → C → A`), no consistent linear ordering exists. Traversing the cycle assigns conflicting positions to the same node (A would be at position 0 and also position 0 + len(A) + len(B) + len(C)).
    -   **Superbubble detection assumes DAGs**: The superbubble algorithm requires acyclic structure. Starting detection from nodes in cycles will produce incorrect results or infinite loops.
    -   **Pre-filtering is efficient**: Detecting all SCCs once (O(V+E)) is cheaper than checking for cycles during every superbubble detection attempt.
-   **How It's Used**:
    -   Nodes with `node_in_scc[v] == true` are **skipped** as superbubble start points in Step 4
    -   `node_component_number` provides **topological ordering** used in tip detection (Step 2) and to ensure superbubbles only chain forward across components (Step 4)
-   **Biological Context**: SCCs often represent tandem repeats, segmental duplications, or assembly errors—regions that genuinely lack a single canonical traversal order.

#### Step 2: Tip Folding (`chainTips`)

-   **Purpose**: Identify and chain "tip" components—dead-end branches that have no forward or backward continuation beyond their component.
-   **Algorithm**:
    -   **Forward tip detection**: Traverse components in reverse topological order. A component is a forward tip if all outgoing edges stay within tip components.
    -   **Backward tip detection**: Traverse components in forward order. A component is a backward tip if all incoming edges originate from tip components.
    -   **Chaining**: Use union-find to merge nodes within tip components via their edges.
-   **Output**:
    -   Union-find structure (`parent`, `rank`) with tip nodes merged
    -   `ignorable_tip[v]`: Boolean flag marking nodes that are part of tips (used to exclude them from subsequent superbubble detection)
-   **Rationale**: Tips represent dead-end branches (e.g., sequence assembly artifacts, singleton variants) that don't participate in variation bubbles. Folding them early prevents them from interfering with superbubble detection.

#### Step 3: Simple Cycle Folding (`chainCycles`)

-   **Purpose**: Detect and collapse simple 2-cycles (bidirectional edges between two nodes).
-   **Algorithm**: For each node `v`:
    -   Find its unique forward neighbor `u_fw` (ignoring tips and self-loops)
    -   Find its unique backward neighbor `u_bw`
    -   If `u_fw == u_bw` and both are unique (not multiple neighbors), then `v` and `u_fw` form a 2-cycle
    -   Union `v` and `u_fw` in the union-find structure and mark `v` as an ignorable tip
-   **Output**: Updated union-find structure and `ignorable_tip` flags
-   **Rationale**: Simple cycles represent bidirectional structural variants or assembly artifacts. Collapsing them simplifies the graph for superbubble detection.

#### Step 4: Superbubble Detection and Chaining (`chainSuperbubbles`)

-   **Purpose**: Identify superbubbles (maximal DAG regions with single entry/exit) and chain nodes within them.
-   **Algorithm**: For each node `s` (excluding nodes in SCCs and ignorable tips):
    1. Run superbubble detection from `s` using the algorithm from [Onodera et al. 2013](https://arxiv.org/pdf/1307.7925) (modified to skip ignorable tips)
    2. If a superbubble `(s, t)` is found and `t` is in a strictly downstream component (to avoid intra-SCC bubbles):
        - Union `s` and `t` in the union-find structure
        - Traverse all nodes reachable from `s` without passing through `t`
        - Union all traversed nodes with `s`, forming a single chain
-   **Superbubble Definition**: A pair of nodes `(s, t)` where:
    -   `s` has out-degree > 0
    -   All paths from `s` reconverge at `t` before reaching any other node
    -   No cycles involve `s` or nodes between `s` and `t`
-   **Output**: Updated union-find structure with superbubble nodes merged into chains
-   **Rationale**: Superbubbles represent variation regions (SNPs, indels, structural variants) where multiple alleles diverge and reconverge. Chaining them enables colinear alignment within variation contexts.

#### Step 5: Chain ID Assignment

-   **Algorithm**:
    1. Path-compress the union-find structure to get stable representative nodes
    2. Collect all unique representatives
    3. Map each representative to a sequential chain ID (0, 1, 2, ...)
    4. Assign each node its chain ID based on its representative
-   **Output**: `node_chain_id[v]` for all nodes

#### Step 6: Linear Coordinate Assignment (`assignLinearPositions`)

-   **Purpose**: Assign approximate linear positions within each chain to enable colinear chaining during mapping.
-   **Algorithm**: For each chain:
    1. Select a deterministic seed node (lowest component number, breaking ties by node ID)
    2. Initialize seed position to a large offset (e.g., `INT32_MAX`) to allow negative coordinates
    3. BFS/DFS propagation:
        - From seed, traverse out-edges: position of successor `u` = position of `v` + length of `v`
        - From seed, traverse in-edges: position of predecessor `u` = position of `v` - length of `u`
    4. Continue until all nodes in the chain have positions assigned
-   **Output**: `node_linear_positions[v]` (int64_t, may be negative)
-   **Properties**:
    -   Nodes along the same path have monotonically increasing positions
    -   Nodes in parallel paths (bubble arms) may have overlapping or interleaved positions
    -   Distance between positions approximates sequence length along edges
-   **Rationale**: Linear coordinates enable O(n log n) colinear chaining via sorting, rather than O(n²) all-pairs distance computation.

#### Final Output

The `AlnGraph` is augmented with per-node metadata:
-   `chain_id`: Chain identifier (0 to num_chains-1). All nodes receive a chain ID (singleton chains for nodes not in bubbles/tips).
-   `linear_coord`: Approximate linear position within the chain (int64_t, may be negative).

This metadata is persisted in the `GraphStore` and queried during the mapping pipeline to accelerate alignment:
-   **Seed clustering**: Group seed hits by `chain_id` in O(n) time
-   **Colinear chaining**: Sort seeds by `linear_coord` and apply gap-cost scoring in O(n log n) time

### Squigglization

This stage converts the base-space sequences of the `AlnGraph` nodes into synthetic signal-space samples.

-   A sliding window of size `k` (from the pore model) moves across each node's sequence.
-   For each k-mer in the window, the pore model provides an expected squiggle value.
-   The result is a sequence of raw signal samples for each node in the `AlnGraph`.

### Quantization

The raw signal samples are processed down two parallel paths for two different purposes: seeding and alignment scoring.

1.  **Fuzzy Quantization**: This process bins the raw signal samples into a smaller set of discrete values. Its goal is to make the subsequent seeding process robust to minor signal variations, allowing similar-looking signals to produce the same seed. This is a lossy conversion designed to improve sensitivity.
2.  **Alignment Quantization**: This process converts the raw signal (typically floats) into the format that will be stored in the `SignalStore` (e.g., `int16_t`). The goal is memory efficiency, while retaining as much precision as possible for the final alignment scoring step.

### Seed Generation & Indexing

This stage builds the `SeedStore`, a hash map for fast seed lookups during mapping.

**Implementation**: `piru/include/index/seed_builder.hpp`, `piru/src/index/seed_builder.cpp`

```cpp
HashSeedStore buildSeedStore(
    const std::vector<FuzzyQuantizedSignal>& signals,
    const SeedExtractor& extractor,
    const SeedBuildConfig& config = SeedBuildConfig{});
```

#### Algorithm Steps

**Step 1: Seed Extraction**

For each node in the graph:
1. Pass the fuzzy-quantized signal to the `SeedExtractor` interface
2. The extractor applies a sliding window and produces a list of seeds
3. Each seed contains: `{hash, position, span}`

**Current extractors**:
- **K-mer extractor** (`KmerSeedExtractor`): Slides a window of `k` consecutive quantized tokens, hashes the window content, and emits seeds at stride intervals
  - Config: `{.backend = "kmer", .k = 10, .stride = 1, .qbits = 4}`
  - Output: Dense seeds every `stride` positions

**Future extractors** (easily added via interface):
- **Minimizer extractor**: Selects local minima within a sliding window to reduce seed density while maintaining coverage
  - Config: `{.backend = "minimizer", .k = 15, .window = 10, .qbits = 4}`
  - Benefit: ~10x fewer seeds, faster mapping, same sensitivity

**Step 2: Hash Table Population**

```cpp
for (node_id in 0..signals.size()) {
    seeds = extractor.extract(signals[node_id])
    for (seed in seeds) {
        store.insert(seed.hash, SeedHit{node_id, seed.position})
    }
}
```

Seeds with the same hash value (collisions) are stored in a list:
```cpp
SeedStore: hash → [(node_id, offset), (node_id, offset), ...]
```

**Step 3: Frequency Statistics**

Compute statistics for seed filtering and scoring:

```cpp
// Max frequency: maximum number of hits for any single hash
max_hash_frequency = max(|hits| for each hash in store)

// Frequency distribution
frequencies = [|hits| for each hash in store]
sort(frequencies)
```

**Step 4: Frequency Filtering (GraphAligner-style)**

To avoid over-represented seeds (repetitive regions, homopolymers), compute a frequency threshold:

```cpp
SeedBuildConfig config{.keep_least_frequent_fraction = 0.1};

// Keep only seeds with frequency ≤ 10th percentile
pos = frequencies.size() * 0.1
threshold = frequencies[pos] + 1
```

**Rationale**: Seeds that occur very frequently (e.g., in repetitive regions) provide little mapping information and slow down alignment. By filtering to the least frequent fraction, we keep informative seeds while discarding noise.

**During mapping**:
```cpp
auto hits = store.lookup(seed_hash);
if (hits && hits->size() < store.frequency_threshold()) {
    // This seed is rare enough to use
    process_seed_hits(*hits);
}
```

#### Extensibility: Switching Seed Types

The `buildSeedStore` function is **backend-agnostic**—it doesn't know whether seeds come from k-mers, minimizers, or any future strategy. To switch seed types:

**Current (k-mer seeds)**:
```cpp
auto extractor = make_seed_extractor(
    SeedExtractorConfig{.backend = "kmer", .k = 10, .stride = 1});
auto store = buildSeedStore(signals, *extractor);
```

**Future (minimizer seeds)**:
```cpp
auto extractor = make_seed_extractor(
    SeedExtractorConfig{.backend = "minimizer", .k = 15, .window = 10});
auto store = buildSeedStore(signals, *extractor);  // NO CHANGE!
```

Adding a new seed type requires:
1. Implement `MinimizerSeedExtractor : public SeedExtractor`
2. Update `make_seed_extractor` factory to recognize `"minimizer"` backend
3. Done! All existing code (hash table, frequency filtering, indexing) works unchanged.

**Design benefit**: Separation of concerns between seed extraction (pluggable strategy) and seed storage (fixed hash table).

#### Output

The `HashSeedStore` provides:
- **Lookup**: `lookup(hash) → vector<SeedHit>*` (returns `nullptr` if hash not found)
- **Frequency stats**: `max_hash_frequency()`, `frequency_threshold()`
- **Size**: `size()` returns unique hash count

**Storage backend**:
```cpp
std::unordered_map<uint64_t, std::vector<SeedHit>> store_;
```

**Alternative backends** (future work):
- Inverted index (CSR format) for memory efficiency
- On-disk hash table for very large graphs
- Filtered hash table (discard seeds above threshold instead of storing flags)

## Index Components (Output)

The final output is a directory containing three main components:

1.  **`GraphStore`**: Stores the topology of the `AlnGraph`. This includes:
    -   Directional nodes and their sequences (the original, un-squigglized bases)
    -   Edges connecting nodes
    -   Per-node chain metadata: `chain_id` and `linear_coord` from pseudo-linearization

    It provides the framework for navigating the graph during alignment and enables efficient seed clustering via chain ID lookups.
2.  **`SignalStore`**: A companion to the `GraphStore`, this stores the **alignment quantized** signal samples for every node's sequence. During mapping, this provides the reference signal to align against.
3.  **`SeedStore`**: The hash table generated from the fuzzy quantized samples. It maps seed hashes to a list of one or more occurrences in the graph, enabling the mapper to quickly find candidate alignment locations.
