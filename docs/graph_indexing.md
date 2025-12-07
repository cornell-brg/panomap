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

**Transformation Strategy: Node Splitting by Successor Context**

The key challenge is that a single VG node may have multiple successors, each providing different k-1 context. We use a **node splitting strategy** to handle this:

-   **Step 1: Bidirectional → Directional Split**: Each bidirectional node is split into two directional nodes (forward and reverse complement), just like DBG transformation.

-   **Step 2: Split by Successor Context**: For each directional node `N` with multiple outgoing edges:
    -   For each successor `S₁, S₂, ..., Sₙ`:
        1. Fetch the first `k-1` bases from successor `Sᵢ`
        2. Create a new node `N_Sᵢ` with sequence = `N.sequence + Sᵢ.sequence[0:k-1]`
        3. Add edge `N_Sᵢ → Sᵢ` with overlap = `k-1`
    -   The original node `N` is replaced by `n` variants: `N_S₁, N_S₂, ..., N_Sₙ`

-   **Step 3: Predecessor Re-wiring**: Update all incoming edges to `N`:
    -   If predecessor `P` had edge `P → N`, it now has edges to all variants: `P → N_S₁, P → N_S₂, ..., P → N_Sₙ`
    -   This preserves all original paths while ensuring each path has correct k-1 context

**Example**:

```
Original VG (no overlaps):
    A (GATTACA) ──→ B (CCC)
                └──→ C (GGG)

After VG→DBG Expansion (k=5, so k-1=4):
    A_B (GATTACA + CCCC) ──→ B (CCC)    [overlap=4]
    A_C (GATTACA + GGGG) ──→ C (GGG)    [overlap=4]

Now squigglization of A_B's last base 'A':
    Sliding window: [ACCC] - has full context from B
```

**Handling Complex Topologies**:

-   **Nodes with single successor**: No splitting needed, just append k-1 bases
-   **Nodes with no successors (tips)**: No expansion needed (or pad with sentinel characters if desired)
-   **Nodes with in-degree > 1**: Each incoming path gets duplicated across all successor splits
-   **Self-loops**: Node gets split to include k-1 bases from itself

**Graph Size Impact**:

-   Worst case: A node with `m` predecessors and `n` successors creates `n` new nodes, each receiving `m` incoming edges
-   Graph size increases, but enables **clean, context-complete squigglization** without runtime edge traversal
-   Trade-off: More nodes vs. simpler squigglization algorithm

**Output**: A DBG-like `AlnGraph` where all nodes are directional and have sequences that include `k-1` bases of overlap from their successors, suitable for squigglization.

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

-   A sliding window is applied to the **fuzzy quantized** samples of each node.
-   The content of the window is hashed to create a "seed." The hashing method is reused from existing signal processing utilities (`piru/signal/`).
-   This seed is stored as a key in a hash table, with its value being a list of locations (`[node_id, offset]`) where the seed occurred in the graph.

## Index Components (Output)

The final output is a directory containing three main components:

1.  **`GraphStore`**: Stores the topology of the `AlnGraph`. This includes:
    -   Directional nodes and their sequences (the original, un-squigglized bases)
    -   Edges connecting nodes
    -   Per-node chain metadata: `chain_id` and `linear_coord` from pseudo-linearization

    It provides the framework for navigating the graph during alignment and enables efficient seed clustering via chain ID lookups.
2.  **`SignalStore`**: A companion to the `GraphStore`, this stores the **alignment quantized** signal samples for every node's sequence. During mapping, this provides the reference signal to align against.
3.  **`SeedStore`**: The hash table generated from the fuzzy quantized samples. It maps seed hashes to a list of one or more occurrences in the graph, enabling the mapper to quickly find candidate alignment locations.
