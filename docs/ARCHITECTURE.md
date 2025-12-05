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
│   (future)   │     │
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
│   (future)   │     │
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
     │   (future)    │
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

**Current status:** Graph loading and model validation are implemented. Graph transformation, pseudo-linearization, squigglization, and index builder are planned.

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
   RawRead         ┌──────────┐
 (raw signal)      │  Index   │
      │            │  Loader  │
      ▼            └──────────┘
 ┌──────────────┐       │
 │Event Detect  │       ├─────── GraphStore
 │ (future)     │       ├─────── SignalStore
 └──────────────┘       └─────── SeedStore
      │
      ▼
 ┌──────────────┐
 │ Normalize    │
 │ (future)     │
 └──────────────┘
      │
      ▼
 Cleaned Signal
      │
      ├────────────────┐
      │                │
      ▼                ▼
  Fuzzy Quant     Alignment Quant
  (for seeding)   (for scoring)
      │                │
      ▼                │
 ┌──────────────┐      │
 │Seed Extract  │      │
 │& Lookup      │      │
 │ (SeedStore)  │      │
 └──────────────┘      │
      │                │
      ▼                │
  Seed Hits            │
 (node_id+pos)         │
      │                │
      ├────────────────┘
      │
      ▼
 ┌──────────────┐
 │   Cluster    │
 │ (by chain ID)│
 └──────────────┘
      │
      ▼
 ┌──────────────┐
 │    Chain     │
 │ (colinear)   │
 └──────────────┘
      │
      ▼
 ┌───────────────┐
 │Align & Score  │
 │(vs SignalStore│
 │ using aln.    │
 │ quant signal) │
 └───────────────┘
      │
      ▼
 AlignmentResult
      │
      ▼
 ┌──────────────┐
 │ ResultWriter │
 │(GAF/GAM/JSON)│
 └──────────────┘
```

**Current status:** Read parsing is implemented. Index loading, signal preprocessing, alignment core, and result writing are planned.

**Mapping algorithm (planned):**
1. **Signal preprocessing**:
   - **Event detection**: Segment raw signal into events or use raw samples directly
   - **Normalize**: Scale/shift signal to standard range
   - Result: Cleaned signal ready for downstream processing
2. **Parallel quantization paths** (cleaned signal branches into two paths):
   - **Path A - Fuzzy quantization**: Apply fuzzy binning (e.g., rawhash2-style) to encourage collision of similar signals for seed discovery
   - **Path B - Alignment quantization**: Convert to format compatible with SignalStore (e.g., float32→int16 if SignalStore uses int16)
3. **Seed extraction & lookup** (using fuzzy quantized signal):
   - Extract signal seeds (k-mer or minimizer) from fuzzy quantized signal
   - Lookup in SeedStore → candidate seed hits (node_id, offset)
4. **Seed clustering** (using GraphStore chain IDs):
   - Group seed hits by chain ID (from pseudo-linearization)
   - Enables O(n) clustering by chain rather than O(n²) all-pairs comparisons
5. **Colinear chaining** (within each cluster):
   - Use linear coordinates from GraphStore for efficient colinear chaining
   - Apply gap costs and scoring to find high-scoring anchor chains
   - Seeds on nodes without chain IDs are handled separately or filtered
6. **Alignment scoring** (using alignment quantized signal):
   - Compare alignment quantized query signal against SignalStore reference signals along chains
   - Produces final alignment coordinates, scores, and statistics

**Note on two quantizations:**
- **Fuzzy quantization** (step 2A): Intentional binning to make similar signals collide for seed discovery. Lossy by design.
- **Alignment quantization** (step 2B): Format conversion to match SignalStore representation (e.g., float32→int16). Lossless or near-lossless for scoring accuracy.
- **Storage quantization** (SignalStore, done during indexing): Compression for memory efficiency on disk/in-memory.

**Note on pseudo-linearization:**
- Chain IDs and linear coordinates (assigned during indexing) enable efficient seed clustering (step 4) and colinear chaining (step 5)

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
- Target type: `AlnGraph` (planned) - directional, cleaned graph ready for squigglization, with chain IDs and linear coordinates.
- Interface: `GraphStore` (planned) - abstract interface for graph topology navigation and chain queries.
- Factory: `make_graph_store` (planned) - creates backend based on representation choice.
- Backends (planned):
  - Adjacency list (default) - simple, fast, cache-friendly.
  - CSR (Compressed Sparse Row) - optimized for static graphs, excellent cache locality.
  - Future: libhandlegraph-compatible backend, GPU/distributed representations.
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
- Interface: `SignalStore` (planned) - abstract interface for storing and retrieving per-node expected signal profiles.
- Factory: `make_signal_store` (planned) - creates backend based on quantization choice.
- **Storage Quantization** (goal: memory efficiency while preserving alignment accuracy):
  - Backends (planned):
    - `Float32SignalStore` - full precision (baseline).
    - `Int16SignalStore` - 16-bit integer quantization.
    - `Int8SignalStore` - 8-bit quantization for aggressive compression.
    - Custom fixed-point or adaptive quantization schemes.
  - Note: This is NOT the same as seed quantization; this is about storage format.
- Usage: During alignment, query signal is compared against reference signals from SignalStore for scoring.

### Seed Storage (Index)
- Interface: `SeedStore` (planned) - abstract interface for seed lookup during mapping.
- Factory: `make_seed_store` (planned) - creates backend based on seeding strategy and quantization choice.
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
  - Implement graph transformation (ImportedGraph → AlnGraph: directional, cleaned)
  - Implement pseudo-linearization (AlnGraph → AlnGraph with chain IDs + linear coordinates):
    - Superbubble detection algorithm
    - Chain assignment and linear coordinate generation
    - Handle nodes that don't belong to chains
  - Implement GraphStore interface + backends (adjacency list, CSR) with chain ID and coordinate storage
  - Implement SignalStore interface + backends (float32, int16, int8, custom quantization)
  - Implement squigglization (AlnGraph + pore model → SignalStore)
  - Implement SeedStore interface + backends:
    - Seeding strategies (k-mer, minimizer)
    - Seed quantization backends (rawhash2-style, other binning schemes)
- **Map**:
  - Implement index loader (GraphStore + SignalStore + SeedStore)
  - Implement alignment core (AlignSingle):
    - Seed lookup via SeedStore
    - Seed clustering by chain ID
    - Colinear chaining within chains using linear coordinates
    - Alignment scoring against SignalStore
  - Integrate ResultWriter
- **Graph I/O**: GBZ loader; richer GFA support; alignment graph definition and export for debugging.
- **Reads**: Investigate POD5 backend.
- **Results**: Richer tags; SAM output support.
