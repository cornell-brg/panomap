# PIRU Architecture

Piru is an adaptive sampling classifier that maps raw nanopore signal directly
to pangenome graphs. It is a classifier first, mapper second -- the goal is
keep/reject decisions on reads based on target regions.

## Subcommands

| Command | Purpose |
|---------|---------|
| `piru index` | Build index from GFA/VG graph + pore model -> .pirx |
| `piru map` | Map reads against index, optionally classify with ROI |
| `piru annotate` | Project BED intervals onto graph -> ROI node set (.pira) |

## Pipeline Overview

### Index
```
GFA/VG + pore model
       |
       v
  GraphLoader -> Transform -> Linearize -> Squigglize -> Build
       |
       v
     .pirx (GraphStore + SeedStore + LinearizationCoords)
```

### Map
```
BLOW5 + .pirx
    |
    v
  Mapper: DSP -> seeds -> lookup -> NodeAnchors
                                        |
                                        v
  Chainer: [expand -> merge -> DP chain] -> ChainResult
                                                |
                                                v
  Mapper: [ROI classify (optional)] -> ResultFormatter -> PAF/GAF
```

### Annotate
```
BED + GFA -> project onto ref path -> expand to haplotypes (mini DP)
          -> walk intervals -> union node sets -> .pira
```

## Responsibility Split

**Mapper** owns: DSP (event detection, normalization, quantization, seeding),
seed lookup, ROI classification, result formatting.

**Chainer** owns: everything between NodeAnchors and ChainResult. The chainer
interface is intentionally opaque -- backends decide their own coordinate space
and internal strategy.

```cpp
class Chainer {
  virtual ChainResult chain(const vector<NodeAnchor>& hits) const = 0;
  virtual string name() const = 0;
};
```

Two backends available:
- **PathChainer** (`--chainer path-chain`, default): expand NodeAnchors to
  per-path PathAnchors (linear space), merge overlapping anchors, DP colinear
  chain per path, extract top-K chains. Fast, cache-friendly, but produces
  alias chains across haplotypes.
- **GraphChainer** (`--chainer graph-chain`): chain directly on NodeAnchors in
  graph space. Cross-path DP with per-path binary search for predecessors.
  True haplotype hopping, deduplicates alias chains, ~1.5x slower.

## Key Types

| Type | Space | Role |
|------|-------|------|
| `SeedEntry` | index | Hash table record: node_id + offset |
| `NodeAnchor` | graph | Seed hit with read context: node_id + offset + read_pos + span |
| `PathAnchor` | linear | Projected anchor: path_id + ref_coord (internal to DPChainer) |
| `ChainedAnchor` | output | Survived DP chaining: has score + chain_id |
| `ChainResult` | output | Best chain(s) with scores and anchor lists |

## Module Map

```
src/commands/       CLI entry points (index.cpp, map.cpp, annotate.cpp)
src/mapping/        Mapping pipeline (batch_mapper, dp_chainer, result_formatter)
src/index/          Graph loading, transformation, linearization, indexing
src/signal/         DSP: event pipelines, fuzzy quantizers, seed extractors
src/io/             Read/result/model I/O (BLOW5, PAF, GAF, pore models)
src/concurrency/    Thread pool (TBB backend)
src/util/           Timing, logging, metrics

include/            Headers (mirrors src/ layout)
tests/              Unit tests (doctest + CTest)
```

## Index Format (.pirx)

Single-directory index:
- **GraphStore**: topology, paths, node sequences
- **SeedStore**: hash table mapping signal seeds -> (node_id, offset) entries
- **Linearization coords**: node -> list of (path_id, ref_coord) for anchor expansion

See `docs/index_format.md` for the binary format spec.

## ROI Classification

For adaptive sampling, reads are classified as keep/reject based on target regions:

1. `piru annotate` projects BED intervals onto the pangenome graph, capturing
   all haplotype-equivalent nodes -> .pira file
2. `piru map --roi file.pira --mode enrich` classifies reads using one of two
   mutually exclusive modes:

**`--chain-genome <overlap>`** -- whole-genome chaining, classify by ROI overlap:
- Chain all anchors normally, then check what fraction of the chain overlaps
  ROI nodes in query space
- Overlap >= threshold -> KEEP, else REJECT (flipped for `--mode deplete`)

**`--chain-target <score>`** -- ROI-only chaining, classify by chain score:
- Filter anchors to ROI nodes before chaining (4-5x faster)
- Classify by chain score >= threshold (overlap is meaningless since all
  anchors are ROI nodes by definition)
- Off-target reads produce weak noise chains below threshold -> REJECT

Output tags: `ro:f:` (overlap or score), `rd:A:K/R` (keep/reject decision).
