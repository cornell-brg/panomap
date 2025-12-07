# Graph Transformation (ImportedGraph → AlnGraph)

## Overview

The graph transformation stage converts the input graph (GFA or VG format) into a directional `AlnGraph` with k-1 overlaps suitable for squigglization. This document focuses on the **DBG (de Bruijn Graph) transformation** implemented in `index/transform_dbg.hpp/cpp`.

## Goals

1. **Convert bidirectional to directional**: Split each node into forward and reverse complement nodes
2. **Preserve k-1 overlaps**: Ensure every node has sufficient overlap context for squigglization
3. **Maintain graph topology**: Preserve all paths from the original graph
4. **Track original mapping**: Store metadata for GAF/GAM output (original node IDs, strand orientation)

## The k-mer Context Problem

### Why Overlaps Matter

Squigglization converts base sequences to expected Nanopore signals by sliding a k-mer window across the sequence and looking up each k-mer in a pore model:

```
Sequence: GATTACA
k=5 window:
  [GATTA] -> lookup -> signal values
   [ATTAC] -> lookup -> signal values
    [TTACA] -> lookup -> signal values
```

At node boundaries, the window needs `k-1` bases from the successor node:

```
Node 1: GATTACA
Node 2: TAG

Without overlap, can't squigglize the last 'A' of Node 1:
  [TTAC?] -> missing context!

With k-1 overlap, Node 1 becomes GATTACA + TAG[0:k-1]:
  [TTACA] -> complete context!
```

### DBG vs Pore Model k-mers

- **DBG k-mer size (`graph_k`)**: Node overlap in the input graph (e.g., 63 for Bifrost)
  - DBG nodes overlap by `graph_k - 1` bases (e.g., 62bp)
- **Pore model k-mer size (`pore_k`)**: Window size for signal lookup (e.g., 5 for R9.4)
  - Only need `pore_k - 1` overlap for squigglization (e.g., 4bp)

Since `pore_k < graph_k`, we **trim** the redundant overlap to save space.

## DBG Transformation Algorithm

### Input
- `ImportedGraph`: Bidirectional graph loaded from GFA
- `graph_k`: DBG k-mer size (e.g., 63)
- `pore_k`: Pore model k-mer size (e.g., 5)

### Output
- `AlnGraph`: Directional graph with trimmed sequences and k-1 overlaps

### Algorithm Steps

#### 1. Validation
```cpp
if (pore_k > graph_k) {
    throw std::invalid_argument("pore_k cannot exceed graph_k");
}
```

**Rationale**: The pore model cannot require more overlap than the DBG provides.

#### 2. Calculate Trimming Amount
```cpp
const std::size_t k_delta = graph_k - pore_k;  // e.g., 63 - 5 = 58
```

**Rationale**: Nodes have `graph_k - 1` overlap, but we only need `pore_k - 1` overlap, so trim the last `k_delta` bases.

#### 3. Node Splitting and Trimming

For each node `N` in the imported graph:

```cpp
// Check node is long enough after trimming
if (N.sequence.size() <= k_delta) {
    LOG_WARN("Node too short, skipping");
    continue;
}

// Trim last k_delta bases
trimmed_seq = N.sequence.substr(0, N.sequence.size() - k_delta);

// Create forward node
AlnNode fwd;
fwd.label = N.id;           // Original GFA node name
fwd.original_id = N.id;     // For GAF/GAM mapping
fwd.is_reverse = false;
fwd.sequence = trimmed_seq;

// Create reverse node
AlnNode rev;
rev.label = N.id;
rev.original_id = N.id;
rev.is_reverse = true;
rev.sequence = revcomp(trimmed_seq);  // Reverse complement of trimmed sequence
```

**Key points**:
- Each bidirectional node becomes two directional nodes (F/R)
- Sequences are trimmed **before** reverse complementing
- Original node ID is preserved for result mapping

#### 4. Edge Rewiring and Overlap Adjustment

For each edge `E` in the imported graph:

```cpp
// Map edge to directional nodes based on orientation
from_id = E.from_reverse ? reverse_node : forward_node;
to_id = E.to_reverse ? reverse_node : forward_node;

// Adjust overlap: original overlap (graph_k - 1) minus trimming amount
std::size_t overlap = E.overlap_bases.value_or(graph_k - 1);
if (overlap >= k_delta) {
    overlap -= k_delta;  // e.g., 62 - 58 = 4
} else {
    overlap = 0;  // Safety
}

graph.addEdge({from_id, to_id, overlap});
```

**Key points**:
- Orientation flags (`from_reverse`, `to_reverse`) determine which directional node to use
- Edge overlaps are reduced by `k_delta` to match trimmed sequences
- Default overlap is `graph_k - 1` if not specified in GFA

### Example

**Input GFA (graph_k=63)**:
```
S   1   GATTACA[...63 bases total]
S   2   TAGCCGC[...63 bases total]
L   1   +   2   +   62M
```

**After transformation (pore_k=5, k_delta=58)**:
```
AlnGraph nodes:
  Node 0 (fwd): GATTACA[...5 bases total]  (original_id=1, is_reverse=false)
  Node 1 (rev): TGTAATC[...5 bases total]  (original_id=1, is_reverse=true)
  Node 2 (fwd): TAGCCGC[...5 bases total]  (original_id=2, is_reverse=false)
  Node 3 (rev): GCGGCTA[...5 bases total]  (original_id=2, is_reverse=true)

Edges:
  0 → 2, overlap=4  (62 - 58)
```

## Implementation Details

### Reverse Complement

```cpp
std::string revcomp(const std::string& seq) {
    // A↔T, C↔G, handles both upper and lowercase
    // Unknown bases (N, etc.) → 'N'
}
```

**Rationale**: Reverse complement is applied to the **trimmed** sequence to maintain correct overlap semantics.

### Edge Orientation Mapping

GFA edges specify orientation with `+` (forward) or `-` (reverse):

```
L   1   +   2   -   62M
     ^       ^
     |       |
from_reverse  to_reverse
```

Mapping to directional nodes:
- `from_reverse=false, to_reverse=false`: forward_1 → forward_2
- `from_reverse=false, to_reverse=true`: forward_1 → reverse_2
- `from_reverse=true, to_reverse=false`: reverse_1 → forward_2
- `from_reverse=true, to_reverse=true`: reverse_1 → reverse_2

### Original Node Mapping

Each `AlnNode` stores:
- `original_id`: The original GFA node name/ID
- `label`: Human-readable label (usually same as original_id)
- `is_reverse`: Boolean indicating strand

This enables GAF/GAM output:
```
// During alignment output
const AlnNode& node = graph.node(aln_id);
gaf_writer.write({
    .node_name = node.label,              // "s1"
    .strand = node.is_reverse ? '-' : '+',
    .position = alignment_position,
});
```

## Edge Cases

### 1. Nodes Too Short After Trimming

If `node.sequence.size() <= k_delta`, the node would be empty or negative length after trimming.

**Handling**: Skip the node with a warning.

**Example**: If `graph_k=63`, `pore_k=5` (k_delta=58), any node shorter than 59bp is skipped.

### 2. Missing Overlap Values

Some GFA files don't specify overlap lengths in `L` records.

**Handling**: Default to `graph_k - 1` (the standard DBG overlap).

### 3. Self-Loops

A node with an edge to itself.

**Handling**: The edge is preserved in both forward and reverse orientations as needed.

## Validation

The function performs two levels of validation:

1. **Input validation**: `pore_k <= graph_k` (throws exception)
2. **Runtime warnings**: Nodes too short after trimming (logs warning, skips node)

## Performance Considerations

- **Time complexity**: O(V + E) where V = nodes, E = edges
- **Space complexity**: O(2V + E) (doubled nodes, same edge count)
- **Memory usage**: Trimming reduces per-node memory by `k_delta` bases

For a graph with 1M nodes and `k_delta=58`, trimming saves ~58MB of sequence data.

## Relationship to VG Transformation

The DBG transformation assumes overlaps already exist. For variation graphs (VG), which lack inherent overlaps, a different transformation is needed (see DEV006):

- **DBG**: Preserve existing overlaps, trim redundant bases
- **VG**: Generate overlaps via node splitting by successor context (Strategy A)

Both produce the same output type (`AlnGraph` with k-1 overlaps), enabling a unified squigglization pipeline.

---

## Squigglization and Quantization

After graph transformation and pseudo-linearization, the next step is to generate expected Nanopore signals from the graph sequences. This process is called **squigglization**.

### Overview

**Implementation**: `piru/include/index/squigglize.hpp`, `piru/src/index/squigglize.cpp`

```cpp
SquiggleResult squigglizeAndQuantize(
    const AlnGraph& graph,
    const io::KmerModel& model,
    const signal::FuzzyQuantizer& fuzzy_quantizer,
    const signal::AlignmentQuantizer& alignment_quantizer);
```

**Input**:
- `AlnGraph`: Transformed graph with k-1 overlaps
- `KmerModel`: Pore model mapping k-mers to expected signal levels (e.g., R9.4, R10.4)
- `FuzzyQuantizer`: For seed extraction (e.g., RH2)
- `AlignmentQuantizer`: For storage (e.g., int16 or passthrough)

**Output**:
- `SquiggleResult`: Two parallel vectors indexed by node ID
  - `fuzzy_signals`: Fuzzy-quantized signals for seed extraction
  - `alignment_signals`: Alignment-quantized signals for storage/mapping

### Algorithm Steps

#### 1. K-mer Sliding Window (First Pass)

For each node, slide a k-mer window across the sequence and look up the expected signal level:

```cpp
for (node_id in 0..graph.nodeCount()):
    sequence = graph.node(node_id).sequence

    if (sequence.size() < k):
        continue  // Skip short nodes

    raw_signals[node_id] = []
    for (i in 0..(sequence.size() - k + 1)):
        kmer = sequence[i:i+k]
        mean_level = model.lookup(kmer)  // pA signal level
        raw_signals[node_id].push_back(mean_level)
        sum += mean_level
        count++
```

**Example (k=5, R9.4 model)**:
```
Sequence: GATTACA (length=7)
K-mers:   GATTA (→ 65.2 pA), ATTAC (→ 72.1 pA), TTACA (→ 68.5 pA)
Raw signal: [65.2, 72.1, 68.5]
```

**Note**: Missing k-mers (containing N or invalid bases) default to 0.0 and log a warning.

#### 2. Global Normalization (Two-Pass Algorithm)

**Pass 1: Compute global mean**
```cpp
global_mean = sum / count
```

**Pass 2: Compute variance** (numerically stable)
```cpp
variance = 0.0
for (node_id in 0..graph.nodeCount()):
    for (val in raw_signals[node_id]):
        diff = val - global_mean
        variance += diff * diff

variance /= count
global_std = sqrt(variance)
```

**Why two-pass?**
- The one-pass formula `variance = E[X²] - E[X]²` suffers from catastrophic cancellation when values are large
- Two-pass algorithm `variance = Σ(x - μ)² / N` is numerically stable

**Apply z-score normalization with outlier clipping**:
```cpp
for each raw signal value:
    z_score = (value - global_mean) / global_std

    // Clip outliers to [-3, 3] range
    if (z_score < -3.0): z_score = -3.0
    if (z_score > 3.0): z_score = 3.0

    normalized_value = z_score
```

**Example**:
```
Raw signals across all nodes: [65.2, 72.1, 68.5, 70.0, 66.8, ...]
Global mean: 68.5 pA
Global std: 3.2 pA

Normalization:
  65.2 → (65.2 - 68.5) / 3.2 = -1.03
  72.1 → (72.1 - 68.5) / 3.2 = 1.13
  68.5 → (68.5 - 68.5) / 3.2 = 0.00
```

#### 3. Quantization

Apply two types of quantization for different purposes:

```cpp
for (node_id in 0..graph.nodeCount()):
    normalized_signal = normalize(raw_signals[node_id], global_mean, global_std)

    // For seed extraction (fuzzy)
    fuzzy_signals[node_id] = fuzzy_quantizer.quantize(normalized_signal)

    // For alignment/storage
    alignment_signals[node_id] = alignment_quantizer.quantize(normalized_signal)
```

**Fuzzy quantization (RH2 example)**:
- Maps normalized signals to reduced alphabet for efficient seed extraction
- Example: RH2 reduces 256-level signals to 2-bit (4-level) alphabet

**Alignment quantization (int16 example)**:
- Converts float32 to int16 for compact storage
- Passthrough option: keep as float32 for maximum precision

### Why Global Normalization?

**Q: Why normalize across all nodes instead of per-node?**

**A: Comparability across the graph**

During read alignment, we compare read signals to graph signals:
```
Read signal (normalized):     [-0.5, 1.2, 0.3, -1.1, ...]
Graph signal (normalized):    [-0.4, 1.3, 0.2, -1.0, ...]
                               ↑ Can directly compute edit distance
```

If each node was normalized independently:
```
Node A (mean=65, std=2): [0.0, 1.5, -1.0]  // z-score relative to Node A
Node B (mean=80, std=5): [0.0, 0.5, -0.5]  // z-score relative to Node B
                          ↑ NOT comparable!
```

Global normalization ensures all signals are on the **same scale**, enabling direct comparison during alignment.

### Edge Cases

**1. Empty nodes** (sequence.size() < k):
- Skip squigglization (no k-mers can be extracted)
- Result vectors contain empty signals for these nodes

**2. All values identical** (stddev ≈ 0):
- Normalize all values to 0.0
- Avoids division by zero

**3. Missing k-mers in model**:
- Default to 0.0 signal level
- Log warning: "Missing k-mer in model at node X pos Y"

### Performance

- **Time complexity**: O(Σ sequence_lengths) - linear in total sequence length
- **Space complexity**: O(Σ (sequence_length - k + 1)) - one float per k-mer position
- **Memory optimization**: Use int16 alignment quantizer to reduce storage by 50%

**Example**: 1M nodes, avg 100bp each, k=5
- Total k-mers: ~96M
- Raw signals: 96M × 4 bytes = 384 MB
- Int16 quantized: 96M × 2 bytes = 192 MB (50% savings)

### Integration with Pipeline

```cpp
// After transformation and pseudo-linearization
AlnGraph graph = transformDbg(imported, graph_k, pore_k);
SccResult scc = computeScc(graph);
TipFoldingResult tips = chainTips(graph, scc);
chainCycles(graph, tips);
SuperbubbleResult sb = chainSuperbubbles(graph, scc, tips);
auto chain_ids = assignChainIds(sb.uf);
auto positions = assignLinearPositions(graph, chain_ids, scc);

// Load pore model
auto model = io::loadKmerModel("r9.4_450bps");

// Create quantizers
auto fuzzy_quantizer = signal::createFuzzyQuantizer("rh2");
auto alignment_quantizer = signal::createAlignmentQuantizer("int16");

// Squigglize
auto squiggle_result = squigglizeAndQuantize(graph, model, fuzzy_quantizer, alignment_quantizer);

// Extract seeds from fuzzy signals
auto seeds = extractSeeds(squiggle_result.fuzzy_signals);

// Store alignment signals in index
index.storeSignals(squiggle_result.alignment_signals);
```

---

## References

- Implementation: `piru/include/index/transform_dbg.hpp`, `piru/src/index/transform_dbg.cpp`
- Squigglization: `piru/include/index/squigglize.hpp`, `piru/src/index/squigglize.cpp`
- Legacy reference: `PIRU-workspace/PIRU/src/graph/alignment_graph.cpp` (lines 541-585)
- Design doc: `piru/docs/graph_indexing.md` (lines 237-244: Squigglization)
- Architecture: `piru/docs/ARCHITECTURE.md` (lines 48-119)
