# Results I/O Overview

This doc summarizes the current alignment result structures and writers for developers.

## Data model

- `AlignmentResult` (`include/io/results/result.hpp`):
  - Basic fields: query name/length/start/end, strand, target path/length/start/end, matches, alignment_block_length, MAPQ.
  - Optional payloads: query sequence/quality.
  - Path/mapping detail: `mappings` vector of `{node_id, offset, is_reverse, edits}`, where edits carry `{from_length, to_length, sequence}`.
  - Optional tags: raw strings ready for GAF (`tp:A:P`, etc.).

## Writers

- Interface: `ResultWriter` (`include/io/results/result_writer.hpp`) with `write(const AlignmentResult&)`.
- Factory: `make_result_writer(path)` (`include/io/results/result_writer_factory.hpp`) selects by extension.

### GAF writer

- Location: `include/io/results/gaf_writer.hpp`, `src/io/results/gaf_writer.cpp`.
- Emits TSV GAF lines (no libvgio dependency). Uses core fields + optional tags.
- Chooses via `.gaf` extension.

### GAM writer (vg protobuf)

- Location: `include/io/results/gam_writer.hpp`, `src/io/results/gam_writer.cpp`.
- Uses libvgio `VGAlignmentEmitter` (format `GAM`) to write vg `Alignment` protobufs.
- Requires libvgio (`PIRU_HAS_LIBVGIO`); fails loudly otherwise.
- Conversion: `to_vg_alignment` (`include/src/io/results/alignment_conversion`) builds a `vg::Alignment` from `AlignmentResult` (sequence/quality, mappings/edits).
- Chooses via `.gam` extension.

### JSON writer (vg JSON)

- Location: `include/io/results/json_writer.hpp`, `src/io/results/json_writer.cpp`.
- Uses libvgio `VGAlignmentEmitter` (format `JSON`) to emit vg `Alignment` records as JSONL.
- Same conversion path as GAM; requires libvgio.
- Chooses via `.json` extension.

## Usage sketch

```cpp
#include "io/results/result_writer_factory.hpp"
#include "io/results/result.hpp"

auto writer = piru::io::make_result_writer("output.gaf"); // or .gam/.json
piru::io::AlignmentResult r;
// populate r.query_name, target_path, mappings, edits, etc.
writer->write(r);
```

## Notes / Roadmap

- GAM/JSON require `PIRU_ENABLE_LIBVGIO` and dependencies (protobuf, htslib, jansson); factory returns `nullptr` or logs errors if unavailable.
- GAF writer is libvgio-free and safe for all builds.
- Future: richer optional tags, batching helpers, and integration hooks in mapping pipeline.
