# VG Graph Transformation (Path-Guided with Expansion)

## Overview

The VG (Variation Graph) transformation converts variation graphs into directional `AlnGraph`s with k-1 overlaps suitable for squigglization. Unlike DBG graphs which have inherent overlaps, VG graphs require generating overlaps through a **path-guided approach with local expansion fallback**.

**Implementation**: `piru/include/index/path_guided_transform.hpp`, `piru/src/index/path_guided_transform.cpp`

## The VG Challenge

### Why VG Graphs Need Transformation

Variation graphs represent genomic variation but lack the k-1 overlaps required for squigglization:

```
VG Graph (no overlaps):
  Node 1: GATTACA
  Node 2: CCGG
  Node 3: TTAGC

During squigglization with k=5:
  Node 1 end: ...TACA?  ← Missing 1 base of context!
```

**Solution**: Use embedded haplotype paths to guide context addition, and expand uncovered nodes by traversing successors.

### DBG vs VG Transformation

| Aspect | DBG | VG |
|--------|-----|-----|
| **Input overlaps** | Present (k-1 bases) | Absent |
| **Approach** | Trim redundant overlap | Generate overlap via paths + expansion |
| **Path information** | Not used | Critical for accuracy |
| **Node duplication** | 2× (forward/reverse) | 2× to 24× (depends on path contexts) |
| **Graph size** | Smaller | Larger (due to context variants) |

## High-Level Algorithm

The transformation follows three stages:

```
Stage 1: ImportedGraph (VG, bidirectional) → AlnGraph (directional, no context)
         ↓
Stage 2: Walk embedded paths, add k-1 context from successors
         ↓
Stage 3: Expand uncovered nodes (not on any path) via successor traversal
```

## Stage 1: Node Splitting (Bidirectional → Directional)

### Goal
Convert each bidirectional VG node into forward and reverse complement variants.

### Algorithm

```cpp
AlnGraph graph;
NodeMapping node_mapping;  // original_id → (forward_aln_id, reverse_aln_id)

for each node in imported_graph:
    // Create forward variant
    AlnNode fwd;
    fwd.label = node.id + "_F";
    fwd.original_id = node.id;
    fwd.is_reverse = false;
    fwd.sequence = node.sequence;  // No context yet
    size_t fwd_id = graph.addNode(fwd);

    // Create reverse variant
    AlnNode rev;
    rev.label = node.id + "_R";
    rev.original_id = node.id;
    rev.is_reverse = true;
    rev.sequence = revcomp(node.sequence);  // No context yet
    size_t rev_id = graph.addNode(rev);

    node_mapping[node.id] = {fwd_id, rev_id};
```

### Edge Transformation

Transform bidirectional edges to directional edges:

```cpp
for each edge in imported_graph:
    // Map to directional nodes based on orientation flags
    from_aln_id = edge.from_reverse ? node_mapping[edge.from].second  // reverse
                                    : node_mapping[edge.from].first;  // forward

    to_aln_id = edge.to_reverse ? node_mapping[edge.to].second
                                : node_mapping[edge.to].first;

    graph.addEdge({from_aln_id, to_aln_id, overlap_bases: 0});
```

**Result**: Directional graph structure with no k-1 context yet.

## Stage 2: Path-Guided Context Addition

### Goal
Use embedded haplotype paths to add k-1 context to nodes, creating variants when contexts differ.

### Algorithm

```cpp
for each path in imported_graph.paths:
    AlnPath aln_path;
    aln_path.name = path.name;

    for i in 0..path.steps.size():
        step = path.steps[i];
        current_node_id = step.segment_id;
        is_reverse = step.is_reverse;

        // Get base AlnGraph node (from Stage 1)
        base_aln_id = is_reverse ? node_mapping[current_node_id].second
                                 : node_mapping[current_node_id].first;

        // Determine successor context
        string context = "";
        size_t successor_aln_id = 0;

        if (i + 1 < path.steps.size()):
            successor_step = path.steps[i + 1];
            successor_aln_id = get_aln_id(successor_step);
            successor_node = graph.node(successor_aln_id);
            context = successor_node.sequence.substr(0, k_minus_1);

        // Create variant key: (original_id, is_reverse, context)
        variant_key = current_node_id + "|" + (is_reverse ? "R" : "F") + "|" + context;

        if (variant_map contains variant_key):
            // Reuse existing node with same context
            node_with_context_id = variant_map[variant_key];
        else:
            // Create new variant with this context
            AlnNode new_variant;
            new_variant.original_id = current_node_id;
            new_variant.is_reverse = is_reverse;
            new_variant.sequence = base_sequence + context;  // Add k-1 context!

            node_with_context_id = graph.addNode(new_variant);
            variant_map[variant_key] = node_with_context_id;

            // Wire to successor
            if (successor_aln_id != 0):
                graph.addEdge({node_with_context_id, successor_aln_id, k_minus_1});

        aln_path.steps.push_back({node_with_context_id, is_reverse});
        covered_nodes.insert(current_node_id);

    graph.addPath(aln_path);
```

### Node Deduplication

**Key insight**: Nodes with the same successor context share the same variant.

**Example**:
```
Path 1: Node A → Node B → Node C
Path 2: Node A → Node B → Node D

At Node B:
- Path 1: successor = C, context = C[0:k-1]
- Path 2: successor = D, context = D[0:k-1]

If context(C) != context(D):
  → Create 2 variants of Node B
Otherwise:
  → Reuse the same Node B variant
```

### Coverage Tracking

While walking paths, track which original nodes were visited:

```cpp
covered_nodes = set();
for each path:
    for each step in path:
        covered_nodes.insert(step.segment_id);

uncovered_nodes = all_nodes - covered_nodes;
```

## Stage 3: Local Expansion for Uncovered Nodes

### Goal
Add k-1 context to nodes not covered by any path by traversing their successors.

### Challenge

Uncovered nodes have no path guidance, so we must explore successor paths to collect k-1 bases. However, deep branching can cause **exponential explosion** if we explore all paths.

**Solution**: Greedy collection - for each immediate successor, follow the *first available path* to collect k-1 bases.

### Algorithm: Context Collection

```cpp
vector<ContextInfo> collectKMinus1Contexts(AlnGraph graph,
                                           size_t start_node_id,
                                           size_t depth,
                                           unordered_set<size_t> visited) {
    vector<ContextInfo> contexts;
    successors = graph.outgoing(start_node_id);

    if (successors.empty()):
        return [{"", 0}];  // Tip node, no context

    // For each immediate successor
    for successor_id in successors:
        if (successor_id in visited):
            continue;  // Cycle detection

        // Greedily collect k-1 bases along first available path
        collected_context = "";
        current_node = successor_id;
        path_visited = visited + {successor_id};

        while (collected_context.size() < depth):
            current_sequence = graph.node(current_node).sequence;

            // Take as many bases as needed
            bases_needed = depth - collected_context.size();
            bases_to_take = min(bases_needed, current_sequence.size());
            collected_context += current_sequence.substr(0, bases_to_take);

            if (collected_context.size() >= depth):
                break;

            // Move to first available successor
            next_successors = graph.outgoing(current_node);
            found_next = false;
            for next_id in next_successors:
                if (next_id not in path_visited):
                    current_node = next_id;
                    path_visited.insert(next_id);
                    found_next = true;
                    break;

            if (!found_next):
                break;  // Cycle or dead end

        contexts.push_back({collected_context, successor_id});

    return contexts;
}
```

**Key points**:
- **Greedy**: Takes first available path at each branch
- **Cycle-safe**: Tracks visited nodes to avoid infinite loops
- **Depth-limited**: Stops after collecting k-1 bases
- **Per-successor**: Creates one context per immediate successor (handles branching)

### Algorithm: Node Expansion

```cpp
void expandUncoveredNodes(AlnGraph graph,
                          set<string> uncovered_node_ids,
                          NodeMapping node_mapping,
                          size_t k_minus_1) {
    for uncovered_id in uncovered_node_ids:
        // Process both forward and reverse variants
        for is_reverse in {false, true}:
            base_aln_id = is_reverse ? node_mapping[uncovered_id].second
                                     : node_mapping[uncovered_id].first;
            base_node = graph.node(base_aln_id);

            // Collect all k-1 contexts from successors
            visited = {};
            contexts = collectKMinus1Contexts(graph, base_aln_id, k_minus_1, visited);

            // Deduplicate contexts by string
            unique_contexts = {};
            for ctx in contexts:
                if (ctx.context not in unique_contexts):
                    unique_contexts[ctx.context] = ctx;

            // Create variant for each unique context
            for (context_str, context_info) in unique_contexts:
                AlnNode variant;
                variant.original_id = base_node.original_id;
                variant.is_reverse = base_node.is_reverse;
                variant.sequence = base_node.sequence + context_str;  // Add k-1 context!

                variant_id = graph.addNode(variant);

                // Wire predecessors to variant
                for pred_id in graph.incoming(base_aln_id):
                    graph.addEdge({pred_id, variant_id, 0});

                // Wire variant to successor
                if (!context_str.empty() && context_info.successor_aln_id != 0):
                    graph.addEdge({variant_id, context_info.successor_aln_id, k_minus_1});
}
```

### Example: Branching Expansion

```
Uncovered node: Node X (sequence: AAAA)
Successors: Node Y (GGGG), Node Z (TTTT)
k=5 (k-1 = 4)

Expansion:
  Variant 1: AAAA + GGGG = AAAAGGGG  → edges to Node Y
  Variant 2: AAAA + TTTT = AAAATTTT  → edges to Node Z

Result: 2 variants of Node X, each with appropriate context
```

### Edge Cases

**1. Tip nodes (no successors)**:
```cpp
if (successors.empty()):
    // Create variant with empty context (just base sequence)
    variant.sequence = base_node.sequence;
```

**2. Cycles**:
```cpp
// Detected via visited set
if (next_id in path_visited):
    break;  // Stop traversal
```

**3. Short sequences**:
```cpp
// May not collect full k-1 bases if paths are too short
// Result: partial context (better than none)
```

## Complete Pipeline Integration

```cpp
AlnGraph PathGuidedTransform::apply(ImportedGraph imported,
                                     size_t graph_k,  // unused for VG
                                     size_t pore_k) {
    // Stage 1: Node splitting
    auto [graph, node_mapping] = importedGraphToAlnGraph(imported);

    // Stage 2: Path-guided context addition
    unordered_set<string> covered_nodes;
    walkPathsAndAddContext(graph, imported, node_mapping, pore_k, covered_nodes);

    // Stage 3: Expand uncovered nodes
    unordered_set<string> uncovered_nodes;
    for (node in imported.nodes):
        if (node.id not in covered_nodes):
            uncovered_nodes.insert(node.id);

    if (!uncovered_nodes.empty()):
        size_t k_minus_1 = pore_k - 1;
        expandUncoveredNodes(graph, uncovered_nodes, node_mapping, k_minus_1);

    // Update statistics
    stats_.original_node_count = imported.nodes.size();
    stats_.transformed_node_count = graph.nodeCount();
    stats_.uncovered_node_count = uncovered_nodes.size();

    return graph;
}
```

## Graph Size Impact

### Node Duplication

**Best case**: All nodes covered by single path with no branching
- Duplication factor: 2× (forward + reverse only)

**Typical case**: Nodes covered by multiple paths with some different contexts
- Duplication factor: 2× to 6× per original node
- Example: 1000 original nodes → 2000-6000 AlnGraph nodes

**Worst case**: Node on 12 haplotype paths, all with different successors
- Duplication factor: up to 24× (12 paths × 2 directions)
- Rare in practice (most nodes share contexts across haplotypes)

**Uncovered nodes**: Each immediate successor creates one variant
- Duplication factor: 2× to 2×N where N = number of successors

### Memory Considerations

For a VG graph with 100K nodes:
- **Before**: 100K nodes × 50bp avg = 5MB sequence data
- **After (typical 4× duplication)**: 400K nodes × (50bp + 4bp context) = 21.6MB
- **Increase**: ~4.3× memory usage

This is acceptable given that:
1. Enables accurate path-specific squigglization
2. Reduces graph size from initial VG representation
3. Allows unified pipeline with DBG graphs

## Performance Characteristics

### Time Complexity

- **Stage 1 (Node splitting)**: O(V + E) where V = nodes, E = edges
- **Stage 2 (Path walking)**: O(P × L) where P = paths, L = avg path length
- **Stage 3 (Expansion)**: O(U × S × k) where U = uncovered nodes, S = avg successors, k = pore_k
- **Overall**: O(V + E + P×L + U×S×k)

### Space Complexity

- **Node storage**: O(V × D) where D = avg duplication factor (2-6×)
- **Edge storage**: O(E × D)
- **Temporary structures**: O(V) for mappings and visited sets

### Practical Performance

For `drb1.vg` (HLA DRB1 locus):
- Input: ~500 nodes, 12 haplotype paths
- Output: ~2000 nodes (4× duplication)
- Transform time: <100ms
- Memory: <10MB

## Validation and Testing

### Unit Tests

**Test suite**: `piru/tests/test_path_guided_transform.cpp`

1. **Test 1-5 (Stage 2)**: Path-guided transformation
   - Single linear path
   - Shared nodes with same context (reuse)
   - Shared nodes with different contexts (duplicate)
   - No paths (all uncovered)
   - Coverage statistics accuracy

2. **Test 6-9 (Stage 3)**: Local expansion
   - Single successor expansion
   - Multiple successors (branching)
   - Tip nodes (no successors)
   - Moderate branching

**All 9 tests passing** ✅

### Integration Testing

```bash
# Test with drb1.vg
./piru index -g tests/data/graphs/drb1.vg -o drb1.index -k 5

# Verify output graph structure
./piru inspect drb1.index
```

## Comparison with DBG Transformation

| Feature | DBG Transform | VG Transform |
|---------|---------------|--------------|
| **Overlaps** | Present, trimmed | Generated from paths |
| **Path usage** | Not used | Critical guidance |
| **Complexity** | O(V + E) | O(V + E + P×L) |
| **Output size** | 2× input | 2-6× input |
| **Accuracy** | Perfect (overlaps given) | Path-dependent |
| **Use case** | Bifrost, cDBG | VG graphs, GBWTs |

## Future Work (Out of Scope)

### Alternative Fallback Strategies (Stage 4)

Currently not implemented, but designed for:

1. **Ambiguous context**: Average k-mer values from multiple successors
2. **Hole marking**: Mark regions with insufficient context as "holes"
3. **Skip**: Simply skip uncovered nodes with warnings

### Full Branch Exploration

Current implementation uses greedy traversal. Future work could explore **all** successor paths:

```cpp
// Recursive exploration (exponential!)
function collectAllContexts(node, depth, visited):
    if depth == 0:
        return [""]

    all_contexts = []
    for successor in node.successors:
        for sub_context in collectAllContexts(successor, depth-1, visited):
            all_contexts.append(successor[0] + sub_context)

    return all_contexts
```

**Trade-off**: More accurate contexts vs exponential runtime.

## References

- **Implementation**: `piru/include/index/path_guided_transform.hpp`, `piru/src/index/path_guided_transform.cpp`
- **Tests**: `piru/tests/test_path_guided_transform.cpp`
- **Planning**: `plans/DEV006-index-vg-support.md`
- **VG Investigation**: `dev-docs/vg_to_dbg_investigation.md`
- **Graph Indexing**: `piru/docs/graph_indexing.md`
- **Architecture**: `piru/docs/ARCHITECTURE.md`

---

*Last updated: December 2024 - Stages 1-3 complete*
