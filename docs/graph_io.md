# Graph I/O Overview

This document describes the current graph I/O interfaces and implementations in PIRU. The goals are:
- keep parsing/serialization concerns isolated from mapping/indexing logic,
- allow swapping front-end formats (GFA, vg) without touching the rest of the code,
- provide a simple in-memory container for downstream conversion into the alignment graph.

## Interfaces

- `ImportedGraph` (`include/io/graphs/graph.hpp`):
  - Lightweight container for parsed graphs.
  - Nodes: `ImportedGraphNode { id, sequence }`.
  - Edges: `ImportedGraphEdge { from, to, from_reverse, to_reverse, overlap, overlap_bases }`:
    - `overlap` holds the raw field: CIGAR for GFA links, stringified length for vg (default `"0"`).
    - `overlap_bases` is the parsed length when available (counts M/= /X for GFA CIGARs; vg overlap length).
  - Paths: `ImportedPath { name, steps, overlaps }` with steps as `ImportedPathStep { segment_id, is_reverse }`.
  - `flavor` tag: `ImportedGraphFlavor` (unknown/dbg/vg) set by the caller (e.g., index flag) to inform downstream processing.
  - Designed as a staging representation; mapping code will later convert this into an internal alignment graph type.

- `GraphLoader` (`include/io/graphs/graph_loader.hpp`):
  - Interface with `bool load(ImportedGraph&)` and `std::string get_format_name() const`.
  - Implementations hide format-specific details; callers use the factory.

- `make_graph_loader` (`include/io/graphs/graph_loader_factory.hpp`):
  - Picks a loader based on file extension; returns `nullptr` if unsupported.

## Implementations

- GFA loader (`src/io/graphs/gfa_loader.cpp`):
  - Supports GFA1 `S` (segments), `L` (links), and `P` (paths).
  - Paths capture oriented segment lists and optional per-edge overlaps.
  - Minimal validation; malformed lines log errors and continue.
  - Edge overlap: preserves CIGAR string; `overlap_bases` counts M/= /X ops when present.

- vg loader (`src/io/graphs/vg_loader.cpp`):
  - Uses libvgio (`VGio::vgio`) to parse vg protobuf streams.
  - Extracts nodes, edges (orientation via `from_start`, `to_end`), and paths (mapping positions).
  - Deduplicates path names and concatenates steps when a path is split across multiple vg Graph messages.
  - Edge overlap: stores numeric overlap length (string and parsed size); defaults to 0 when unspecified.

## Using the loader

```cpp
#include "io/graphs/graph_loader_factory.hpp"

auto loader = piru::io::make_graph_loader(graph_path);
if (!loader) { /* handle unsupported format */ }
piru::io::ImportedGraph graph;
if (!loader->load(graph)) { /* handle parse failure */ }
// graph.nodes / graph.edges / graph.paths now available for conversion

// Example: iterate elements
for (const auto& node : graph.nodes) {
    // node.id, node.sequence
}
for (const auto& edge : graph.edges) {
    // edge.from, edge.to, edge.from_reverse, edge.to_reverse
    // edge.overlap (raw), edge.overlap_bases (parsed length if available)
}
for (const auto& path : graph.paths) {
    // path.name, path.overlaps (raw per-edge strings if present)
    for (const auto& step : path.steps) {
        // step.segment_id, step.is_reverse
    }
}
```

The `index` command currently uses this flow to load a graph, log basic stats, and will later transform it into the alignment graph for indexing.

## Dependencies

- GFA: no external deps.
- vg: requires libvgio (and its deps: protobuf, htslib, jansson). Enabled via `PIRU_ENABLE_LIBVGIO` CMake option; defines `PIRU_HAS_LIBVGIO` when linked.

## Roadmap

- Add GBZ (GBWT + Graph) support once available; likely via an additional loader behind the same factory.
- Broaden GFA support (walks, containments, optional fields) and stricter validation.
- Enhance vg loader with richer path metadata (offsets, edits) if needed by the alignment graph builder.
