# PIRU Index File Format Specification

**Version**: 1.0
**Status**: Draft
**Last Updated**: 2025-01-07

## Overview

PIRU indexes are stored as a multi-file directory structure following the VG-style approach. Each index consists of three binary files, one for each major component:

- **`.graph`**: Graph topology, chain metadata, and linear coordinates (GraphStore)
- **`.signals`**: Per-node signal profiles for alignment scoring (SignalStore)
- **`.seeds`**: Seed hash table for rapid candidate lookup (SeedStore)

Each file is independently loadable and uses a common header structure for versioning and backend identification.

## Directory Layout

```
index_name/
├── index_name.graph      # GraphStore (topology + chain metadata)
├── index_name.signals    # SignalStore (per-node signal profiles)
└── index_name.seeds      # SeedStore (seed hash table)
```

**Naming conventions:**
- Index directory name: user-specified via `--output` (default: `<graph-file>.piru`)
- Component files: `basename.{graph,signals,seeds}` where `basename` matches the directory name

**Example:**
```bash
$ piru index -o my_genome.piru input.gfa
# Creates: my_genome.piru/my_genome.{graph,signals,seeds}
```

## Common Header Structure

All three files share a common header format for versioning, integrity checking, and backend identification.

### Header Layout

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        8      char[8]    magic              Magic number (file type identifier)
8        4      uint32_t   format_version     Format version (major.minor encoded as major*1000 + minor)
12       4      uint32_t   header_size        Total header size in bytes (for forward compatibility)
16       4      uint32_t   backend_length     Length of backend type string
20       N      char[]     backend_type       Backend type string (e.g., "adjlist", "int16", "hash")
20+N     M      ...        backend_metadata   Backend-specific metadata (optional, backend-defined)
```

**All multi-byte integers are stored in little-endian format.**

### Magic Numbers

Each file type has a unique 8-byte magic number:

| File Type | Magic (ASCII)    | Magic (Hex)               |
|-----------|------------------|---------------------------|
| `.graph`  | `PIRUGRAF`       | `50 49 52 55 47 52 41 46` |
| `.signals`| `PIRUSIGS`       | `50 49 52 55 53 49 47 53` |
| `.seeds`  | `PIRUSEED`       | `50 49 52 55 53 45 45 44` |

### Format Version

Encoded as: `major * 1000 + minor`

Example: Version 1.0 → `1000`
Example: Version 1.2 → `1002`

**Versioning policy:**
- **Major version bump**: Breaking changes (incompatible file structure)
- **Minor version bump**: Backward-compatible additions (new optional fields)

### Backend Type String

Identifies the storage backend used for this component. Enables polymorphic loading via factory pattern.

**Common backend types:**

| Component   | Backend Type | Description                          |
|-------------|--------------|--------------------------------------|
| GraphStore  | `adjlist`    | Adjacency list representation        |
|             | `csr`        | Compressed sparse row (future)       |
| SignalStore | `float32`    | 32-bit floating point signals        |
|             | `int16`      | 16-bit integer quantized signals     |
|             | `int8`       | 8-bit integer quantized (future)     |
| SeedStore   | `hash`       | Hash table (unordered_map)           |
|             | `csr`        | CSR inverted index (future)          |

## GraphStore Format (`.graph` file)

The `.graph` file stores the directional alignment graph topology along with pseudo-linearization metadata (chain IDs and linear coordinates).

### File Structure

```
[Common Header]
[Global Index Metadata]
[Graph Metadata]
[Node Array]
[Edge Array]
[Path Array] (optional, may be empty)
```

### Global Index Metadata

This section stores build-time parameters and provenance information needed by the mapper to ensure compatibility.

```
Offset   Size   Type       Field                  Description
------   ----   ----       -----                  -----------
H        4      uint32_t   piru_version_major     PIRU major version
H+4      4      uint32_t   piru_version_minor     PIRU minor version
H+8      4      uint32_t   piru_version_patch     PIRU patch version
H+12     8      uint64_t   build_timestamp        Unix timestamp (seconds since epoch)
H+20     4      uint32_t   graph_flavor           0=unknown, 1=dbg, 2=vg
H+24     4      uint32_t   graph_k                Graph k-mer size (DBG parameter)
H+28     4      uint32_t   pore_k                 Pore model k-mer size
H+32     4      uint32_t   model_name_length      Length of pore model name string
H+36     N      char[]     model_name             Pore model name (e.g., "r9.4_450bps")
H+36+N   4      uint32_t   fuzzy_quantizer_length Length of fuzzy quantizer backend name
H+40+N   M      char[]     fuzzy_quantizer        Fuzzy quantizer type (e.g., "rh2")
H+40+N+M 4      uint32_t   align_quantizer_length Length of alignment quantizer backend name
H+44+N+M P      char[]     align_quantizer        Alignment quantizer type (e.g., "int16")
H+44+N+M+P 4    uint32_t   source_path_length     Length of source graph path
H+48+N+M+P Q    char[]     source_path            Original graph file path (optional, for provenance)
```

Where `H` = header size from common header.

**Graph flavor values:**
- `0` = Unknown
- `1` = DBG (de Bruijn graph)
- `2` = VG (variation graph)

**Notes:**
- PIRU version allows the mapper to warn about version mismatches
- `model_name` is the canonical model identifier (e.g., "r9.4_450bps", "r10.4_400bps")
- `fuzzy_quantizer` and `align_quantizer` ensure mapper uses compatible quantization
- `source_path` is optional and may be relative or absolute (for debugging/provenance)

### Graph Metadata

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
G        8      uint64_t   node_count         Number of nodes
G+8      8      uint64_t   edge_count         Number of edges
G+16     8      uint64_t   path_count         Number of paths
G+24     4      uint32_t   reserved           Reserved (must be 0)
G+28     4      uint32_t   reserved           Reserved (must be 0)
```

Where `G` = offset after Global Index Metadata section.

### Node Array

Array of `node_count` node records:

```
Per-node record:
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   sequence_length    Length of sequence in bases
4        N      char[]     sequence           DNA sequence (ACGT, no null terminator)
4+N      4      uint32_t   chain_id           Chain ID from pseudo-linearization
8+N      8      int64_t    linear_position    Linear coordinate within chain
16+N     4      uint32_t   label_length       Length of label string
20+N     M      char[]     label              Node label (e.g., "n1_fwd", no null terminator)
```

**Notes:**
- Nodes are stored in order (node 0, node 1, ..., node N-1)
- Node IDs are implicit (array index)
- Sequence is stored as ASCII characters ('A', 'C', 'G', 'T')

### Edge Array

Array of `edge_count` edge records:

```
Per-edge record:
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   from_node          Source node ID
4        4      uint32_t   to_node            Target node ID
8        4      uint32_t   overlap            Overlap length in bases
```

**Notes:**
- Edges are directional (from_node → to_node)
- Overlap indicates how many bases the sequences share at the junction
- No duplicate edges (GraphStore validates uniqueness)

### Path Array (Optional)

Array of `path_count` path records. May be empty if no paths are stored.

```
Per-path record:
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   name_length        Length of path name
4        N      char[]     name               Path name (no null terminator)
4+N      4      uint32_t   step_count         Number of steps in path
8+N      M      ...        steps              Array of step_count step records

Per-step record (within path):
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   node_id            Node ID visited
4        1      uint8_t    is_reverse         0 = forward, 1 = reverse
```

### Backend-Specific Metadata (Adjacency List)

For `backend_type = "adjlist"`, no additional backend metadata is stored in the header. The adjacency structure is implicit from the edge array.

## SignalStore Format (`.signals` file)

The `.signals` file stores per-node quantized signal profiles used for alignment scoring.

### File Structure

```
[Common Header]
[Signal Metadata]
[Signal Array]
```

### Signal Metadata

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
H        8      uint64_t   node_count         Number of nodes (must match .graph)
H+8      4      uint32_t   quantization_bits  Bits per sample (8, 16, or 32)
H+12     4      float      scale              Quantization scale factor
H+16     4      float      offset             Quantization offset
H+20     4      uint32_t   reserved           Reserved (must be 0)
```

**Quantization formula:**
- **Encoding**: `quantized = round((raw_signal - offset) / scale)`
- **Decoding**: `raw_signal = quantized * scale + offset`

### Signal Array

Array of `node_count` signal records:

```
Per-signal record:
Offset   Size   Type         Field              Description
------   ----   ----         -----              -----------
0        4      uint32_t     sample_count       Number of signal samples
4        N      int8/16/32[] samples            Quantized signal samples
```

**Sample types by backend:**
- `int8`: 8-bit signed integers (`int8_t`)
- `int16`: 16-bit signed integers (`int16_t`)
- `float32`: 32-bit floats (`float`)

**Notes:**
- Signals are stored in node order (signal for node 0, node 1, ..., node N-1)
- Sample count may vary per node (depends on sequence length and pore model k)
- For a sequence of length L and pore model k, sample count = L - k + 1

### Backend-Specific Metadata (int16)

For `backend_type = "int16"`:

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      int16_t    min_value          Minimum representable value (-32768)
4        4      int16_t    max_value          Maximum representable value (32767)
```

## SeedStore Format (`.seeds` file)

The `.seeds` file stores the seed hash table mapping seed hashes to lists of (node_id, offset) hit tuples.

### File Structure

```
[Common Header]
[Seed Metadata]
[Hash Entry Array]
[Hit List Array]
```

### Seed Metadata

```
Offset   Size   Type       Field                  Description
------   ----   ----       -----                  -----------
H        8      uint64_t   unique_hash_count      Number of unique seed hashes
H+8      8      uint64_t   total_hit_count        Total number of seed hits across all hashes
H+16     4      uint32_t   max_hash_frequency     Maximum hits for any single hash
H+20     8      uint64_t   frequency_threshold    Frequency filter threshold
H+28     8      double     filter_fraction        Keep least frequent fraction (0.0-1.0)
H+36     4      uint32_t   extractor_name_length  Length of the seed extractor name string
H+40     N      char[]     extractor_name         Name of the seeder (e.g., "kmer", "minimizer")
H+40+N   4      uint32_t   num_extractor_params   Number of key-value parameters for the extractor
```
This is followed by `num_extractor_params` records of the following structure:
```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   key_length         Length of parameter key string
4        N      char[]     key                Parameter key (e.g., "k", "stride", "window_size")
4+N      4      uint32_t   value_length       Length of parameter value string
8+N      M      char[]     value              Parameter value (as a string)
```

### Hash Entry Array

Array of `unique_hash_count` hash entry records:

```
Per-hash entry:
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        8      uint64_t   hash               Seed hash value
8        8      uint64_t   hit_offset         Offset into hit list array
16       4      uint32_t   hit_count          Number of hits for this hash
```

**Notes:**
- Entries are stored in ascending hash order for binary search
- `hit_offset` is byte offset into the Hit List Array section
- Hashes with zero hits (filtered) are not stored

### Hit List Array

Concatenated array of all seed hits, indexed by hash entries:

```
Per-hit record:
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
0        4      uint32_t   node_id            Node containing this seed
4        4      uint32_t   offset             Position of seed within node's signal
```

**Example:**
- Hash `0x123ABC` has 3 hits at `hit_offset=0, hit_count=3`
- The 3 hit records start at byte offset 0 in the Hit List Array
- Next hash starts at offset `0 + 3*8 = 24`

### Backend-Specific Metadata (hash)

For `backend_type = "hash"`, no additional backend metadata is needed beyond the common seed metadata.

## Loading Procedure

### Index Loader Workflow

1. **Open index directory** and verify all three files exist
2. **Read and validate headers** from each file:
   - Check magic numbers match expected values
   - Verify format versions are compatible
   - Extract backend type strings
3. **Read and validate global metadata** from `.graph` file:
   - Extract PIRU version, build timestamp
   - Extract build parameters: `graph_k`, `pore_k`, `model_name`
   - Extract quantizer types: `fuzzy_quantizer`, `align_quantizer`
   - **Warn if PIRU version mismatch** (different major version)
   - **Verify pore model compatibility** with mapper's available models
4. **Validate consistency**:
   - `node_count` must match between `.graph` and `.signals`
   - `graph_k` and `pore_k` must match expected values
   - Quantizer types must be supported by mapper
5. **Instantiate backends** via factories:
   - `make_graph_store(backend_type)` → read `.graph` data
   - `make_signal_store(backend_type)` → read `.signals` data
   - `make_seed_store(backend_type)` → read `.seeds` data
6. **Return index components and metadata** to caller

### Pseudo-code

```cpp
struct IndexMetadata {
    uint32_t piru_version_major;
    uint32_t piru_version_minor;
    uint32_t piru_version_patch;
    uint64_t build_timestamp;
    uint32_t graph_flavor;
    uint32_t graph_k;
    uint32_t pore_k;
    std::string model_name;
    std::string fuzzy_quantizer;
    std::string align_quantizer;
    std::string source_path;
};

struct LoadedIndex {
    IndexMetadata metadata;
    std::unique_ptr<GraphStore> graph;
    std::unique_ptr<SignalStore> signals;
    std::unique_ptr<SeedStore> seeds;
};

LoadedIndex load_index(const std::string& index_dir) {
    // Open files
    auto graph_file = open(index_dir + "/" + basename + ".graph");
    auto signal_file = open(index_dir + "/" + basename + ".signals");
    auto seed_file = open(index_dir + "/" + basename + ".seeds");

    // Read headers
    auto graph_hdr = read_common_header(graph_file, "PIRUGRAF");
    auto signal_hdr = read_common_header(signal_file, "PIRUSIGS");
    auto seed_hdr = read_common_header(seed_file, "PIRUSEED");

    // Read global metadata from .graph file
    auto metadata = read_global_metadata(graph_file);

    // Validate version compatibility
    if (metadata.piru_version_major != CURRENT_VERSION_MAJOR) {
        LOG_WARN("Index built with PIRU v" + to_string(metadata.piru_version_major) +
                 ", current version is v" + to_string(CURRENT_VERSION_MAJOR));
    }

    // Validate consistency
    validate_node_counts(graph_hdr, signal_hdr);
    validate_parameters(metadata, seed_hdr);
    validate_quantizers(metadata.fuzzy_quantizer, metadata.align_quantizer);

    // Load via factories
    auto graph = load_graph_store(graph_file, graph_hdr);
    auto signals = load_signal_store(signal_file, signal_hdr);
    auto seeds = load_seed_store(seed_file, seed_hdr);

    return {std::move(metadata), std::move(graph), std::move(signals), std::move(seeds)};
}
```

## Versioning and Compatibility

### Current Version: 1.0

**Format version 1.0 supports:**
- DBG graphs (directional, with overlaps)
- Adjacency list graph representation
- Int16 and float32 signal quantization
- Hash-based seed storage

### Future Compatibility

**Backward compatibility guarantee:**
- Readers must support all minor versions within the same major version
- Example: A reader for v1.2 must be able to read v1.0 and v1.1 files

**Forward compatibility (within major version):**
- Newer readers can skip unknown fields using `header_size`
- Unknown backend types trigger graceful error (cannot load)

**Version upgrade path:**
- v1.x → v2.x: Requires index rebuild (breaking change)
- v1.0 → v1.x: Automatic (backward compatible)

## Optional Features (Future Extensions)

The following features may be added in future minor versions:

### Checksums (v1.1+)

Add optional CRC32 checksums to each file:

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
H-8      4      uint32_t   header_checksum    CRC32 of header (excluding this field)
EOF-4    4      uint32_t   data_checksum      CRC32 of all data after header
```

### Compression (v1.2+)

Add optional compression for SignalStore:

```
Backend metadata (compressed backends):
Offset   Size   Type       Field                  Description
------   ----   ----       -----                  -----------
0        4      uint32_t   compression_type       0=none, 1=zstd, 2=lz4
4        8      uint64_t   uncompressed_size      Original size before compression
12       8      uint64_t   compressed_size        Size of compressed data blob
```

### Memory-Mapped (mmap) Optimization (v1.3+)

Add alignment padding for mmap-friendly layouts:

```
Offset   Size   Type       Field              Description
------   ----   ----       -----              -----------
H        ?      ...        padding            Align data section to page boundary (4KB)
```

## Implementation Notes

### Endianness

**All multi-byte integers use little-endian byte order.**

- This matches x86-64 native byte order (most common platform)
- Cross-platform support (big-endian systems) is out of scope for v1.0
- Future versions may add endianness flag in header

### String Encoding

- All strings are UTF-8 encoded
- No null terminators (length-prefixed instead)
- Labels and names are case-sensitive

### Error Handling

Loaders must validate:
- Magic numbers match expected values
- Format versions are supported
- Node counts are consistent between files
- File size matches expected size from metadata
- All node/edge IDs reference valid nodes

Invalid indexes must trigger clear error messages, not crashes.

### Performance Considerations

- **Sequential reads**: Files are designed for sequential reading (minimize seeks)
- **Cache-friendly**: Node/signal/seed arrays stored contiguously
- **No inter-file dependencies**: Each file is independently loadable (enables parallel loading)

## References

- **Architecture**: `piru/docs/ARCHITECTURE.md` (lines 330-357)
- **Graph Indexing Pipeline**: `piru/docs/graph_indexing.md` (lines 425-437)
- **VG File Formats**: https://github.com/vgteam/vg/wiki/File-Formats (multi-file design reference)
- **DEV005 Ticket**: `plans/DEV005-index-serialization.md`

## Changelog

### Version 1.0 (2025-01-07)
- Initial format specification
- Support for DBG graphs, int16/float32 signals, hash-based seeds
- Multi-file directory layout
- Common header structure with versioning and backend identification
- Global index metadata section with build parameters and provenance:
  - PIRU version (major, minor, patch)
  - Build timestamp
  - Graph flavor (DBG/VG)
  - Graph and pore k-mer sizes
  - Pore model name (e.g., "r9.4_450bps")
  - Fuzzy and alignment quantizer types
  - Source graph path (optional provenance)
