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

Currently one backend:
- **DPChainer** (Method 1): expand NodeAnchors to PathAnchors (linear space),
  merge overlapping anchors, DP colinear chain, extract top-K chains.

Future backends choose their own strategy:
- **Method 2**: seed-hit-centric chaining in query/graph space (haplotype hopping)
- **Method 3**: ROI-only chaining (speed optimization for adaptive sampling)

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
2. `piru map --roi file.pira --mode enrich` maps reads normally, then checks
   chain-ROI overlap in query space
3. Overlap >= threshold -> KEEP, else REJECT (flipped for `--mode deplete`)

Output tags: `ro:f:` (overlap fraction), `rd:A:K/R` (keep/reject decision).

## Vision

**Current (Method 1):** Map against whole genome, DP chain in linear space,
post-check chain overlap with ROI. 94.1% accuracy on yeast benchmark.

**Next (Method 2):** Seed-hit-centric chaining in query/graph space. Each seed
hit carries coords on all paths; consecutive pairs chain if they share any path
with good gap cost. Eliminates alias chains, enables haplotype hopping.

**Future (Method 3):** ROI-only chaining. Only chain seeds landing on ROI nodes.
Off-target reads with zero ROI seeds are instant rejects. Speed optimization --
may lose "negative signal" (strong off-target mapping = confident reject).
