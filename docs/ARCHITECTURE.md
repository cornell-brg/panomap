# PIRU Architecture

## High-Level Command Flows (Blackbox View)

### Index Command
```
   GFA/vg          pore kmer
   graph            model
     │               │
     └───────┬───────┘
             ▼
      ┌─────────────┐
      │ piru index  │
      └─────────────┘
             │
             ▼
          index
         (on disk)
    ┌─────────────────┐
    │ • GraphStore    │
    │ • SignalStore   │
    │ • SeedStore     │
    └─────────────────┘
```

### Map Command
```
    slow5          index
    reads        (on disk)
      │              │
      └──────┬───────┘
             ▼
      ┌─────────────┐
      │  piru map   │
      └─────────────┘
             │
             ▼
        alignments
      (GAF/GAM/JSON)
```

## Detailed Pipeline Breakdowns

Status legend:
- [x] implemented
- [ ] not implemented / planned

### Index Pipeline (Internal Flow)

```
   GFA/vg          pore kmer
   graph            model
     │               │
     ▼               │
┌──────────────┐     │
│ GraphLoader  │     │
└──────────────┘     │
     │               │
     ▼               │
ImportedGraph        │
(bidirectional       │
 raw graph)          │
     │               │
     ▼               │
┌──────────────┐     │
│  Transform   │     │
│              │     │
└──────────────┘     │
     │               │
     ▼               │
  AlnGraph           │
(directional,        │
 cleaned up)         │
     │               │
     ▼               │
┌──────────────┐     │
│Pseudo-Linear │     │
│              │     │
└──────────────┘     │
     │               │
     ▼               │
  AlnGraph           │
(with chain IDs,     │
 linear coords)      │
     │               │
     └───────┬───────┘
             ▼
     ┌───────────────┐
     │ Squigglize &  │
     │ Index Builder │
     │               │
     └───────────────┘
             │
             │ produces three components:
             │
   ┌─────────┼─────────┐
   │         │         │
   ▼         ▼         ▼
GraphStore SignalStore SeedStore
(topology) (per-node   (hash:
            signals)    seed→
                       node+pos)
   │         │         │
   └─────────┴─────────┘
             │
             ▼
       index
     (on disk)
```

**Current status:** Full indexing pipeline is implemented for both DBG and VG graphs (path-guided transform, pseudo-linearization, squigglization/quantization, seed build, serialization).

**Component status:**
- [x] GraphLoader (GFA/vg)
- [x] Transform (ImportedGraph → AlnGraph; DBG transform + VG path-guided transform)
- [x] Pseudo-Linearize (chain IDs, linear coords)
- [x] Squigglize & Index Builder
- [x] GraphStore / SignalStore / SeedStore backends + serialization

**Pipeline stages:**
1. **GraphLoader**: Parses GFA/vg files into ImportedGraph (bidirectional graph with nodes, edges, paths)
2. **Transform**: Converts ImportedGraph → AlnGraph (directional, cleaned up for alignment)
3. **Pseudo-Linearize**: Detects superbubbles (DAGs) in AlnGraph, chains them, and assigns linear coordinates (chain IDs) to nodes. Some nodes may not belong to chains. Enables efficient seed clustering and colinear evaluation during mapping
4. **Squigglize & Index Builder**: Applies pore model to AlnGraph sequences, builds index components

**Index components:**
- **GraphStore**: Stores graph topology (nodes, edges, paths, chain IDs, linear coordinates) for navigation and clustering during alignment. Uses interface + backend pattern (adjacency list, CSR, etc.) to allow experimentation with different representations
- **SignalStore**: Stores per-node expected signal profiles with configurable storage quantization (float32, int16, custom, etc.) for memory efficiency. Uses interface + backend pattern
- **SeedStore**: Hash table mapping signal seeds to (node_id, offset) tuples. Uses interface + backend pattern for:
  - Seeding strategy (k-mer vs minimizer)
  - Seed quantization (binning/bucketing to encourage collision of similar signals, e.g., rawhash2-style)

### Map Pipeline (Internal Flow)

```
    slow5              index
    reads            (on disk)
      │                 │
      ▼                 │
 ┌──────────────┐       │
 │  Read Parse  │       │
 │(ReadProvider)│       │
 └──────────────┘       │
      │                 │
      ▼                 ▼
   RawRead         ┌──────────────┐
 (raw signal)      │  Index       │
      │            │  Loader      │
      ▼            │ (GraphStore, │
 ┌──────────────┐  │  SeedStore,  │
 │Event Detect  │  │  SignalStore,│
 └──────────────┘  │  LinearCoords│
      │            └──────────────┘
      ▼                 │
 ┌──────────────┐       │
 │ Normalize    │       │
 └──────────────┘       │
      │                 │
      ▼                 │
 Cleaned Signal         │
      │                 │
      ├────────────────┐│
      │                ││
      ▼                ▼▼
  Fuzzy Quant     Alignment Quant
  (for seeding)   (for scoring)
      │                │
      ▼                │
 ┌──────────────┐      │
 │Seed Extract  │      │
 └──────────────┘      │
      │                │
      ▼                │
 ┌──────────────┐      │
 │ Seed Lookup  │◀─────┘ (SeedStore)
 └──────────────┘
      │
      ▼
  Seed Hits
(graph space: node_id, offset)
      │
      ▼
 ┌──────────────┐
 │   Anchor     │◀────── (LinearCoords)
 │  Expansion   │
 └──────────────┘
      │
      ▼
    Anchors
(linear space: path_id, ref_coord)
      │
      ▼
 ┌──────────────┐
 │  Clustering  │
 │  /Chaining   │
 └──────────────┘
      │
      ▼
Selected Anchors/Chains
      │
      └─── future: align → score → results ───┘
```

**Current status (DEV012 - 2025-12-15):** Uniform 4-stage mapping pipeline implemented. Seed lookup, anchor expansion, and clustering/chaining are complete. Alignment scoring and result writing remain for future work.

**Component status:**
- [x] Read parse (ReadProvider)
- [x] Event detection
- [x] Signal normalization
- [x] Fuzzy quantization (seeding)
- [x] Alignment quantization (scoring)
- [x] Seed extraction
- [x] Index Loader
- [x] SeedStore lookup
- [x] Anchor expansion (SuperbubbleExpander, PathWalkExpander)
- [x] Clustering/chaining (FSE, Probe, DPChain)
- [ ] Alignment scoring (future)
- [ ] ResultWriter integration in map path (future)

**Mapping algorithm (implemented DEV012):**

**Stage 1: Linearization** (done during indexing)
   - Superbubble: Assign chain IDs and local linear positions to nodes
   - Path-walk: Assign global positions along reference paths

**Stage 2: Signal preprocessing**
   - **Event detection**: Segment raw signal into events or use raw samples directly
   - **Normalize**: Scale/shift signal to standard range
   - **Parallel quantization**: Cleaned signal branches into two paths:
     - **Fuzzy quantization** (for seeding): Apply fuzzy binning to encourage seed collision
     - **Alignment quantization** (for scoring): Convert to SignalStore format

**Stage 3: Seed extraction & lookup**
   - Extract signal seeds from fuzzy quantized signal
   - Lookup in SeedStore → candidate seed hits (graph space: node_id, offset)

**Stage 4: Anchor expansion** (graph → linear transformation)
   - **SuperbubbleExpander**: 1:1 mapping using chain_id from GraphStore
   - **PathWalkExpander**: 1:N mapping via path occurrence coordinates
   - Result: Anchors in linear space (path_id, ref_coord)

**Stage 5: Clustering/Chaining** (select optimal anchor subset)
   - **FSE/Probe clusterers**: Group by path_id, cluster by diagonal (superbubble pipeline)
   - **DPChain clusterer**: Colinear chaining via dynamic programming (path-walk pipeline)
   - Result: Selected anchors/chains for alignment extension

**Stage 6: Alignment scoring** (future work)
   - Compare alignment quantized query signal against SignalStore reference signals
   - Produces final alignment coordinates, scores, and statistics

**Note on quantization:**
- **Fuzzy quantization** (Stage 2): Intentional binning to make similar signals collide for seed discovery. Lossy by design.
- **Alignment quantization** (Stage 2): Format conversion to match SignalStore representation (e.g., float32→int16). Lossless or near-lossless for scoring accuracy.
- **Storage quantization** (SignalStore, done during indexing): Compression for memory efficiency on disk/in-memory.

**Note on linearization:**
- **Superbubble linearization** (Stage 1): Assigns chain IDs and local linear positions to nodes within superbubbles. Fast O(n) clustering. Best for simple variation graphs.
- **Path-walk linearization** (Stage 1): Assigns global positions along reference paths. Enables haplotype-aware chaining. Best for complex graphs with cycles.
- Both linearization strategies enable the same uniform 4-stage pipeline (Stage 2-5)

---

## Module Reference

### CLI and Commands
- Entry: `src/main.cpp` wires subcommands (`index`, `map`, `eval`, `mt-test`).
- Parsing: `cli/parse` turns argv into `Parsed`.
- Handlers: `src/commands/*.cpp` consume args, log summaries, and call into IO layers.

### Read I/O
- Interface: `ReadProvider` (`include/io/reads/read_provider.hpp`).
- Factory: `make_read_provider` (`include/io/reads/read_provider_factory.hpp`).
- Backends: slow5/blow5 (`Slow5Provider`), POD5 stub.
- Data: `RawRead` carries signal and metadata.

### Model I/O
- Interface: `KmerModel` (`include/io/models/model.hpp`).
- Factory: `load_builtin_model`, `load_model_from_file` (`include/io/models/model_factory.hpp`).
- Data: built-in R9.4/R10.4 models + file loader.
- Usage: Pore model is consumed during **index** construction to convert base sequences → expected signal profiles. It is not needed during **map** since mapping works directly with signal data.

### Graph I/O
- Staging type: `ImportedGraph` (`include/io/graphs/graph.hpp`) with nodes/edges/paths, overlaps, and flavor tag (dbg/vg/unknown).
- Interface: `GraphLoader` (`include/io/graphs/graph_loader.hpp`).
- Factory: `make_graph_loader` (`include/io/graphs/graph_loader_factory.hpp`).
- Backends: GFA parser (S/L/P), vg loader via libvgio (nodes/edges/paths).
- Notes: GBZ/GBWTGraph support is on the roadmap; alignment graph conversion will happen downstream.

### Graph Storage (Index)
- Target type: `AlnGraph` - directional, cleaned graph ready for squigglization, with chain IDs and linear coordinates.
- Interface: `GraphStore` (implemented).
- Backends:
  - `AdjListGraphStore` (current default) - simple adjacency-list wrapper over `AlnGraph`.
  - Future: CSR/other compact layouts, libhandlegraph-compatible backend, GPU/distributed representations.
- **Stored data**:
  - Topology: nodes, edges, paths
  - **Chain IDs**: which chain (if any) each node belongs to (from pseudo-linearization)
  - **Linear coordinates**: position within chain for colinear ordering
- **Usage**:
  - Node lookup, edge traversal, path queries during alignment
  - Chain ID lookup for seed clustering (group seeds by chain)
  - Linear coordinate queries for efficient colinear chaining within chains
  - Separate from SignalStore (which holds signal data)

### Signal Storage (Index)
- Interface: `SignalStore` (implemented).
- Backend: `VectorSignalStore` (current) storing alignment-quantized per-node signals.
- Notes: Alignment quantizer backend is selectable during indexing; signals are stored using that quantized format. Future backends can introduce alternative storage layouts or precisions.

### Seed Storage (Index)
- Interface: `SeedStore` (implemented).
- Backend: `HashSeedStore` (current) built from fuzzy-quantized seeds; serialized to disk.
- **Two orthogonal design choices:**
  1. **Seeding Strategy**: How seeds are extracted
     - k-mer: consecutive k signal samples.
     - Minimizer: sparsified sampling for reduced index size.
  2. **Seed Quantization** (goal: fuzzy matching via intentional collision of similar signals):
     - Rawhash2-style: bin signal values so similar signals collide in hash table.
     - Other binning/bucketing strategies to balance sensitivity vs specificity.
- Backends (planned):
  - `KmerSeedStore` + quantization backend (e.g., `Rawhash2Quantizer`).
  - `MinimizerSeedStore` + quantization backend.
- Usage: During mapping, query signal is quantized using the same seed quantization scheme, seeds are extracted, and looked up in the index to find candidate (node_id, offset) hits.

### Index Serialization (On-Disk Format)
- **Approach**: Multiple files per index (VG-style), following the modular component architecture.
- **Directory Structure**:
  ```
  index_name/
  ├── index_name.graph      # GraphStore (topology)
  ├── index_name.signals    # SignalStore (per-node signal profiles)
  └── index_name.seeds      # SeedStore (seed hash table)
  ```
- **File Format** (per component):
  ```
  [Magic number: 8 bytes]     # e.g., "PIRUGRAF", "PIRUSIGS", "PIRUSEED"
  [Format version: 4 bytes]   # versioning for backward compatibility
  [Backend type: variable]    # e.g., "adjlist", "csr", "float32", "int16"
  [Component header: variable] # metadata (sizes, parameters, checksums)
  [Data blob: variable]       # actual index data
  ```
- **Benefits**:
  - Each component independently loadable and swappable
  - Matches modular architecture (can swap GraphStore backend without rebuilding SignalStore)
  - Proven approach (VG uses separate .xg, .gcsa, .gbwt files)
  - Easy to version and experiment with different backends
- **Future Optimizations**:
  - mmap support for zero-copy loading (plan data layouts accordingly)
  - Compression for SignalStore (zstd/lz4) if size becomes an issue
  - Checksums for integrity verification
- **Loading**: Index loader reads all three files and instantiates appropriate backend classes via factories.

### Result I/O
- Data: `AlignmentResult` (`include/io/results/result.hpp`) with query/target spans, mappings/edits, and tags.
- Interface: `ResultWriter` (`include/io/results/result_writer.hpp`).
- Factory: `make_result_writer` (`include/io/results/result_writer_factory.hpp`) chooses by extension.
- Backends:
  - GAF (pure TSV, no libvgio).
  - GAM (libvgio VGAlignmentEmitter).
  - JSON (libvgio VGAlignmentEmitter JSON).
- Conversion: `to_vg_alignment` (`include/io/results/alignment_conversion.hpp`) builds vg::Alignment for GAM/JSON.
- Note: GAF/GAM/JSON formats will capture alignment statistics, quality metrics, and mapping metadata. Unaligned reads and detailed diagnostics may be added as separate outputs in the future.

### Concurrency & Utilities
- Concurrency: `concurrency/executor` with TBB (optional) or serial fallback.
- Timing: `util/timing` with profiling macros.
- Logging: `util/logging`.
- Signal handling: `util/signal_handlers`.
- Metrics: `util/metrics` for CPU time, real time, and peak RSS.

### Dependency Notes
- Third-party: slow5lib (reads), libvgio (vg graphs + GAM/JSON writers), kmer_models (builtin models).
- Build options: `PIRU_ENABLE_LIBVGIO`, `PIRU_USE_TBB`, `PIRU_FETCH_*` flags in `CMakeLists.txt`.

## Future/Roadmap (selected)
- **Index**:
  - Implement CSR backend for GraphStore.
  - Implement advanced storage backends for SignalStore (custom fixed-point, adaptive) beyond the current vector store.
  - Implement minimizer-based seeding strategy for SeedStore and corresponding storage tweaks.
- **Map**:
  - Implement alignment core (AlignSingle):
    - Seed lookup via SeedStore
    - Seed clustering by chain ID
    - Colinear chaining within chains using linear coordinates
    - Alignment scoring against SignalStore
  - Integrate ResultWriter
- **Graph I/O**: GBZ loader; richer GFA support; alignment graph definition and export for debugging.
- **Reads**: Investigate POD5 backend.
- **Results**: Richer tags; SAM output support.
