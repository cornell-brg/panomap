# Index Pipeline Test Coverage Gaps

This document identifies missing test coverage for the indexing pipeline components. Tests marked with L are critical gaps that should be addressed before production use.

## Priority 1: Critical Gaps (Must Add)

### L Chain ID Assignment (`assignChainIds`)

**Current Coverage**: 0/3 critical paths tested

**Missing Tests**:

1. **Path compression and representative finding**
   ```cpp
   TEST_CASE("Chain IDs: path compression finds stable representatives") {
       // Create union-find with multiple merges
       // Verify find() returns same representative for merged nodes
   }
   ```

2. **Deterministic chain ID assignment**
   ```cpp
   TEST_CASE("Chain IDs are deterministic and sequential") {
       // Create graph with multiple chains
       // Verify representatives sorted before ID assignment
       // Verify chain IDs are 0, 1, 2, ... (no gaps)
       // Run twice, verify same IDs produced
   }
   ```

3. **Multiple chains get different IDs**
   ```cpp
   TEST_CASE("Separate chains get distinct chain IDs") {
       // Create graph with 3 disconnected components
       // Verify each component gets different chain ID
   }
   ```

**Why Critical**: Chain IDs are used for seed clustering during mapping. Non-deterministic or incorrect IDs break colinear chaining.

---

### L Linear Coordinate Assignment (`assignLinearPositions`)

**Current Coverage**: 0/6 critical paths tested

**Missing Tests**:

1. **Component-aware root selection**
   ```cpp
   TEST_CASE("Linear positions: root selected by (component, node_id)") {
       // Create chain with nodes in different components
       // Verify root is node with lowest (component_number, node_id)
       // NOT just lowest node_id
   }
   ```

2. **Forward propagation**
   ```cpp
   TEST_CASE("Linear positions propagate forward correctly") {
       // Linear chain: 0 í 1 í 2
       // Node lengths: [10, 20, 30]
       // Verify: pos[1] = pos[0] + 10, pos[2] = pos[1] + 20
   }
   ```

3. **Backward propagation**
   ```cpp
   TEST_CASE("Linear positions propagate backward correctly") {
       // Linear chain: 0 ê 1 ê 2
       // Verify: pos[0] = pos[1] - len[0], pos[1] = pos[2] - len[1]
   }
   ```

4. **Negative coordinates allowed**
   ```cpp
   TEST_CASE("Linear positions can be negative") {
       // Chain with backward propagation from root
       // Verify some positions < 0
   }
   ```

5. **Separate chains have independent coordinate spaces**
   ```cpp
   TEST_CASE("Multiple chains have independent linear coordinates") {
       // Create 2 chains with different chain IDs
       // Both should start from similar offsets (e.g., INT32_MAX)
       // Positions in chain A don't affect chain B
   }
   ```

6. **DFS traversal order**
   ```cpp
   TEST_CASE("Linear positions use DFS (not BFS) for traversal") {
       // Create diamond: root í {A, B} í exit
       // Verify positions assigned via DFS (depth-first, not breadth-first)
   }
   ```

**Why Critical**: This is the most complex pseudo-linearization algorithm. Legacy code had bugs here. Zero test coverage = zero confidence.

---

### L DBG Transform: Edge Orientations

**Current Coverage**: 1/4 orientation combinations tested (only FíF)

**Missing Tests**:

1. **Reverse-to-reverse edges**
   ```cpp
   TEST_CASE("DBG transform: reverse í reverse edge orientation") {
       // Edge: n1(-) í n2(-)
       // Verify edge maps to: n1_rev í n2_rev
   }
   ```

2. **Forward-to-reverse edges**
   ```cpp
   TEST_CASE("DBG transform: forward í reverse edge orientation") {
       // Edge: n1(+) í n2(-)
       // Verify edge maps to: n1_fwd í n2_rev
   }
   ```

3. **Reverse-to-forward edges**
   ```cpp
   TEST_CASE("DBG transform: reverse í forward edge orientation") {
       // Edge: n1(-) í n2(+)
       // Verify edge maps to: n1_rev í n2_fwd
   }
   ```

4. **Original node metadata preserved**
   ```cpp
   TEST_CASE("DBG transform preserves original node metadata") {
       // Verify original_id, label, is_reverse flags correct
       // Forward node: is_reverse = false
       // Reverse node: is_reverse = true
   }
   ```

**Why Critical**: Reverse complement handling is error-prone. Only testing FíF means 75% of edge cases uncovered.

---

## Priority 2: High-Impact Gaps (Should Add)

### † Squigglization: Real K-mer Lookup

**Current Coverage**: Uses `ConstKmerModel` (constant signal), never exercises real lookup

**Missing Tests**:

1. **K-mer lookup with varying signal values**
   ```cpp
   TEST_CASE("Squigglization looks up k-mers in pore model") {
       // Create model with k=3:
       //   AAA í 50.0
       //   CCC í 70.0
       //   GGG í 60.0
       // Sequence: AAACCCGGG
       // Verify raw signals: [50, 50, ..., 70, 70, ..., 60, 60]
       // Verify normalization produces different values
   }
   ```

2. **Outlier clipping to [-3, 3]**
   ```cpp
   TEST_CASE("Squigglization clips outliers to [-3, 3]") {
       // Create signal with extreme outliers
       // Most values: 100.0, one outlier: 1000.0
       // After normalization, verify outlier clamped to 3.0
   }
   ```

3. **Missing k-mers default to 0.0**
   ```cpp
   TEST_CASE("Squigglization handles missing k-mers") {
       // Sequence with N: ACGNGTA (k=3)
       // K-mers: ACG, CGN (missing), GNG (missing), NGT (missing), GTA
       // Verify missing k-mers í 0.0 signal
   }
   ```

4. **Two-pass variance calculation**
   ```cpp
   TEST_CASE("Squigglization uses stable two-pass variance") {
       // Create signals with large values (e.g., 1e6)
       // Verify variance computed correctly (no catastrophic cancellation)
   }
   ```

**Why High-Impact**: Current test bypasses k-mer lookup entirely. No confidence the core functionality works.

---

### † Seed Builder: Lookup Correctness

**Current Coverage**: Frequency stats tested, but not actual lookups

**Missing Tests**:

1. **Lookup returns correct hits**
   ```cpp
   TEST_CASE("Seed store lookup returns correct (node_id, offset)") {
       // Build store from 2 nodes with known seeds
       // Lookup specific hash
       // Verify returned SeedHits have correct node_id and offset
   }
   ```

2. **Hash collisions handled correctly**
   ```cpp
   TEST_CASE("Seed store handles hash collisions") {
       // Create scenario where 2 different positions produce same hash
       // Verify lookup returns BOTH hits
   }
   ```

3. **Lookup non-existent hash returns nullptr**
   ```cpp
   TEST_CASE("Seed store lookup returns nullptr for missing hash") {
       // Build store with known seeds
       // Lookup hash that doesn't exist
       // Verify lookup() returns nullptr (not empty vector)
   }
   ```

4. **Seeds extracted with k > 1**
   ```cpp
   TEST_CASE("Seed extraction works with multi-token k-mers") {
       // Signal: [1, 2, 3, 4, 5]
       // k=3, stride=1
       // Verify seeds at positions 0, 1, 2 (windows: [1,2,3], [2,3,4], [3,4,5])
   }
   ```

**Why High-Impact**: Hash table is the core data structure for mapping. Need to verify lookups actually work.

---

### † DBG Transform: Edge Cases

**Current Coverage**: Only happy path tested

**Missing Tests**:

1. **Nodes too short after trimming**
   ```cpp
   TEST_CASE("DBG transform skips nodes too short after trimming") {
       // k_graph=63, pore_k=5, k_delta=58
       // Node with sequence length 50 (< 59)
       // Verify node skipped (not in output graph)
   }
   ```

2. **Missing overlap metadata defaults to graph_k - 1**
   ```cpp
   TEST_CASE("DBG transform uses default overlap when not specified") {
       // Edge without overlap_bases field
       // Verify overlap = graph_k - 1
   }
   ```

3. **Self-loops preserved**
   ```cpp
   TEST_CASE("DBG transform preserves self-loops") {
       // Node with edge to itself
       // Verify self-loop mapped to forward and reverse nodes
   }
   ```

4. **Multiple edges between same nodes**
   ```cpp
   TEST_CASE("DBG transform handles multiple edges") {
       // n1 í n2 (overlap 4)
       // n1 í n2 (overlap 5) [different orientation or context]
       // Verify both edges preserved
   }
   ```

---

## Priority 3: Nice-to-Have Tests

### Tip Folding: More Coverage

**Current Coverage**: Only forward tips tested

**Missing Tests**:

1. **Backward tips detected**
   ```cpp
   TEST_CASE("Tip folding detects backward tips") {
       // Create graph with backward tip (no incoming edges from non-tips)
       // Verify marked as ignorable_tip
   }
   ```

2. **Union-find merging within tips**
   ```cpp
   TEST_CASE("Tip folding merges nodes within tip components") {
       // Tip component with 3 nodes connected by edges
       // Verify all 3 merged in union-find
   }
   ```

3. **Mixed tip and non-tip successors**
   ```cpp
   TEST_CASE("Tip folding handles nodes with both tip and non-tip edges") {
       // Node with 2 outgoing edges: one to tip, one to main chain
       // Verify node NOT marked as tip (has non-tip successor)
   }
   ```

---

### Superbubble: Complex Cases

**Current Coverage**: Simple diamond bubble tested

**Missing Tests**:

1. **Nested superbubbles**
   ```cpp
   TEST_CASE("Superbubble detection handles nested bubbles") {
       // Outer bubble: 0 í {1, 2} í 5
       // Inner bubble: 1 í {3, 4} í 2
       // Verify both detected and chained separately
   }
   ```

2. **Long superbubbles**
   ```cpp
   TEST_CASE("Superbubble detection handles long paths") {
       // Bubble with 10+ nodes in each arm
       // Verify all interior nodes merged
   }
   ```

3. **Superbubbles across components**
   ```cpp
   TEST_CASE("Superbubble chaining skips intra-SCC bubbles") {
       // Start and exit in same SCC component
       // Verify bubble NOT chained (would create conflicts)
   }
   ```

---

### SCC Detection: Self-loops

**Current Coverage**: 3-node cycle tested

**Missing Tests**:

1. **Self-loops detected as SCCs**
   ```cpp
   TEST_CASE("SCC detection marks self-loops as in_scc") {
       // Node with edge to itself
       // Verify in_scc[node] == true
   }
   ```

2. **Multiple disjoint cycles**
   ```cpp
   TEST_CASE("SCC detection handles multiple separate cycles") {
       // Two separate 3-cycles
       // Verify each cycle gets different component ID
       // Verify component numbers respect topological order
   }
   ```

---

### AlnGraph Validation: Negative Tests

**Current Coverage**: Implementation exists, but never called in tests

**Missing Tests**:

1. **Validation detects out-of-bounds edges**
   ```cpp
   TEST_CASE("AlnGraph validate detects invalid edge targets") {
       // Create graph with 3 nodes
       // Add edge 0 í 999 (out of bounds)
       // Verify validate() returns false
   }
   ```

2. **Validation detects asymmetric adjacency**
   ```cpp
   TEST_CASE("AlnGraph validate detects missing back-edges") {
       // Add edge 0 í 1 to out_edges[0]
       // Don't add 0 to in_edges[1] (asymmetry)
       // Verify validate() returns false
   }
   ```

3. **Validation detects inconsistent metadata**
   ```cpp
   TEST_CASE("AlnGraph validate detects metadata inconsistency") {
       // Set chain_id but not linear_position
       // Verify validate() returns false
   }
   ```

---

## Coverage Summary

| Component | Current Tests | Missing Critical Tests | Confidence |
|-----------|---------------|------------------------|------------|
| SCC Detection | 1 | 2 self-loop/multi-cycle cases | † Medium |
| Tip Folding | 1 | 3 backward/union-find cases | L Low |
| Cycle Folding | 1 | 0 (adequate) |  Good |
| Superbubble | 1 | 3 nested/long/cross-component | † Medium |
| **Chain ID Assignment** | **0** | **3 determinism/uniqueness** | L **None** |
| **Linear Positions** | **0** | **6 propagation/root-selection** | L **None** |
| DBG Transform | 1 | 7 orientations/edge-cases | L Low |
| Squigglization | 1 | 4 real-lookup/clipping | L Low |
| Seed Builder | 1 | 4 lookup-correctness | L Low |
| AlnGraph Validation | 0 | 3 negative tests | † Medium |
| **Total** | **7** | **~35** | L **~20%** |

---

## Recommended Testing Workflow

1. **Add Priority 1 tests first** (Chain IDs, Linear Positions, DBG orientations)
   - These are critical algorithms with zero or inadequate coverage
   - Should be added before any production use

2. **Add Priority 2 tests** (Squigglization real lookup, Seed lookups)
   - Current tests use mocks that bypass core functionality
   - High risk of silent failures in production

3. **Add Priority 3 tests** (Edge cases, negative tests)
   - Increases robustness and confidence
   - Can be deferred if time-constrained

4. **Add end-to-end integration test** (Stage 8 of DEV004)
   - Run full pipeline on real GFA file
   - Verify all components integrate correctly

---

## Test Quality Guidelines

**Good tests should**:
-  Use **realistic inputs** (not just trivial cases)
-  Verify **multiple properties** (not just "doesn't crash")
-  Test **edge cases** (empty inputs, boundary conditions)
-  Include **negative tests** (verify failures caught correctly)
-  Use **deterministic data** (avoid random values)

**Avoid**:
- L Mock objects that bypass core logic (like `ConstKmerModel`)
- L Tests that only check output size (not correctness)
- L Trivial inputs (1-2 nodes, single edge)
- L Checking only "happy path" (need failure cases too)

---

## References

- Current tests: `piru/tests/test_index_pipeline.cpp`
- Pipeline documentation: `piru/docs/graph_indexing.md`
- Pseudo-linearization details: `piru/docs/pseudo_linearization.md`
- DEV004 plan: `plans/DEV004-index-dbg-pipeline.md`
