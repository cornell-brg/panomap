# PANOMAP Architecture

Piru maps raw nanopore signal directly to pangenome graphs without basecalling.

## Subcommands

| Command | Purpose |
|---------|---------|
| `panomap index` | Build index from GFA graph + pore model -> .pirx |
| `panomap map` | Map reads against index -> GAF/PAF |

## Pipeline Overview

### Index
```
GFA + pore model
       |
       v
  GraphLoader -> Transform -> Path-walk bucket indexer
       |
       v
     .pirx (GraphStore + SeedStore + LinearizationCoords)
```

Per-path build flow:
```text
for each embedded path
  walk path
    -> squigglize (pore model k-mer lookup)
    -> diff filter (remove near-duplicate events)
    -> tokenize (rh2: adaptive quantization, landmark: peak rise/fall)
    -> extract seeds (kmer or minimizer sliding window)
    -> map seed position back to (node_id, offset)
    -> append hit to bucket[hash & mask]

after all paths
  for each bucket
    sort by (hash, node_id, offset)
    dedup exact duplicates
    build bucket hash table
```

### Map
```
BLOW5 + .pirx
    |
    v
  [per chunk]:
    event detect + normalize -> diff filter -> tokenize -> seed extract
    -> lookup -> accumulate hits -> chain -> mark survivors -> prune
    -> [next chunk or exit]
    |
    v
  final chain -> mapping decision -> GAF output
```

## Signal Processing Pipeline

Two tokenizer backends convert normalized signal to discrete tokens for seeding:

**RH2 (default):** Adaptive quantization of event values into 4-bit bins.
One token per event. Seeds = k consecutive tokens hashed (default k=6).

**Landmark:** Detects valleys (landmarks) in the signal, extracts peaks between
consecutive valleys, encodes each peak by log-quantized rise/fall amplitudes
(2 bits each = 4-bit token). Seeds = k consecutive peak tokens hashed (k=4).
Produces ~3x fewer seeds, ~4x fewer hits, with equivalent accuracy.

The diff filter runs before both tokenizers, removing events within a threshold
of the previous emitted event. For landmark, the prominence filter provides
additional noise rejection.

Position tracking: `original_positions` maps token indices back to event-space
coordinates through both diff filter compression and tokenizer compression.
Seed extractors compute event-space span from these positions.

## Chaining

**Mapper** owns: signal processing, seed lookup, mapping decision, output.

**Chainer** owns: everything between NodeAnchors and ChainResult.

```cpp
class Chainer {
  virtual ChainResult chain(const vector<NodeAnchor>& hits) const = 0;
  virtual string name() const = 0;
};
```

Three backends:
- **PathChainer** (`--chainer path-chain`, default): expand NodeAnchors to
  per-path PathAnchors (linear space), merge overlapping anchors, DP colinear
  chain per path, extract top-K chains.
- **SortChainer** (`--chainer sort-chain`): 1D sort-based chaining with
  pre-computed canonical coordinates. O(anchors) time.
- **PanChainer** (`--chainer pan-chain`): 1D-banded cross-path colinear
  chaining. Combines SortChainer scaling with PathChainer accuracy.

## Key Types

| Type | Space | Role |
|------|-------|------|
| `SeedEntry` | index | Bucketed hash-table hit record: node_id + offset |
| `NodeAnchor` | graph | Seed hit with read context: node_id + offset + read_pos + span |
| `ChainedAnchor` | output | Survived DP chaining: ref_coord + node info |
| `ChainResult` | output | Best chain(s) with scores and anchor lists |

## Module Map

```
src/commands/       CLI entry points (index.cpp, map.cpp)
src/mapping/        Mapping pipeline (batch_mapper, chainers, seed_merger)
src/index/          Graph loading, transformation, linearization, indexing
src/signal/         DSP: event pipeline, diff filter, tokenizers, seed extractors
src/io/             Read/result/model I/O (BLOW5, GAF, PAF, pore models)
src/concurrency/    Thread pool (TBB backend)
src/util/           Timing, logging, metrics, tracing

include/            Headers (mirrors src/ layout)
tests/              Unit tests (doctest + CTest)
```

## Index Format (.pirx)

Single-file binary index:
- **GraphStore**: topology, paths, node sequences
- **SeedStore**: bucket hash table mapping signal seeds -> (node_id, offset)
- **Linearization coords**: node -> list of (path_id, ref_coord) for anchor expansion
- **1D coordinates**: optional canonical node positions for SortChainer/PanChainer

The SeedStore is bucket-finalized and read-only:
- low hash bits choose a bucket
- `keys[]` are sorted inside each bucket for binary search
- exact duplicates removed on `(hash, node_id, offset)` during build
- singleton hashes occupy one entry slot directly
- multi-hit hashes point to a span in the bucket's flat `entries[]` array

See `docs/index_format.md` for the binary format spec.

## GAF Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `pn:Z:` | string | Path name (when graph walk is in col 6) |
| `tp:A:` | char | P=primary, S=secondary |
| `cs:i:` | int | Chain score |
| `an:i:` | int | Anchor count |
| `se:f:` | float | Score per event span (score / query_span) |
| `ad:f:` | float | Anchor density (anchors / ref_span) |
| `ck:i:` | int | Chunks processed (primary only) |
| `dt:f:` | float | Processing time in seconds (primary only) |
