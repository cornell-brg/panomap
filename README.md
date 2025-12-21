<p align="center">
  <img src="docs/assets/logo.png" alt="PIRU" width="200">
</p>

---

<p align="center">
  <em>Real-time squiggle-to-graph mapping for Nanopore adaptive sampling</em>
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

## What is PIRU?

PIRU maps raw Nanopore signal directly to pangenome graphs without basecalling. By working at the squiggle level, it enables real-time decision-making during adaptive sampling experiments. The tool supports multiple graph formats (GFA, vg), SLOW5/BLOW5 signals, and pore chemistries (R9.4, R10.4).

## Quick Start

**Requirements:** CMake 3.16+, C++20 compiler, zlib, TBB (or use `-DPIRU_FETCH_TBB=ON`)

```bash
# Clone and build
git clone --recursive https://github.com/xxx/piru.git
cmake -S piru -B piru/build && cmake --build piru/build

# Create index from a graph
./piru/build/piru index -o my_index reference.gfa

# Map reads
./piru/build/piru map -i my_index reads.blow5
```

**Optional dependencies:**
- VG support: `libjansson-dev`, `libhts-dev` (or `-DPIRU_ENABLE_LIBVGIO=OFF` to skip)
- Tests: `ctest --test-dir piru/build`

## Usage

<details>
<summary><code>piru --help</code></summary>

```
piru 0.0.1 (experimental scaffold)

Usage: piru [--help] [--version] <command> [options]
Subcommands:
  index    Build index (dbg default, --graph=vg for vg).
  map      Map reads against an index (stub).
  eval     Evaluate theoretical enrichment/depletion (stub).
  mt-test  Spawn parallel sleep tasks to test concurrency.

Run 'piru <command> --help' for command-specific options.
```
</details>

<details>
<summary><code>piru index --help</code></summary>

```
Usage: piru index [options] <graph-file>

Required arguments:
  <graph-file>      Graph file to index

Options:
  -h, --help            Show help
  -g, --graph           Graph type: dbg (default) or vg
  -k, --graph-k         DBG k-mer size (default: auto-detect from overlap)
  -m, --model           Pore model (builtin name: r9.4/r10.4 or model file path)
  -o, --output          Output index directory (default: <graph-file>.piru)
  -t, --threads         Worker threads
  -p, --profile         Emit timing profile (tree)

Seed Generation Options:
  --seed-k              Seed k-mer size (default: 6)
  --seed-stride         Seed stride (default: 1)
  --seed-filter         Keep least frequent seed fraction (default: 0.9)

Alignment Quantization Options:
  --aq-backend          Backend: int16 (default), int8, passthrough
  --aq-scale            Manual scale override (expert)
```
</details>

<details>
<summary><code>piru map --help</code></summary>

```
Usage: piru map [options] (--index <dir> | --graph <file> --model <model>) <reads-path>

Required arguments:
  <reads-path>       Input slow5/blow5 file or directory containing reads

Options:
  -h, --help            Show help
  -i, --index           Path to pre-built index directory (XOR with --graph)
  -g, --graph           Graph file for in-memory indexing (XOR with --index)
  -m, --model           Pore model (required with --graph; r9.4/r10.4 or file path)
  -t, --threads         Worker threads (-1 = auto)
  -p, --profile         Emit timing profile (tree)

In-Memory Indexing Options (with --graph):
  --linearizer          Linearizer backend: superbubble (default) or path-walk
  --graph-type          Graph type: dbg (default) or vg
  --graph-k             DBG k-mer size (default: auto-detect from overlap)

Mapping Options:
  --max-seed-freq       Maximum seed frequency for lookup (default: use index threshold)
  --clusterer           Clusterer backend: fse (default), probe, dp-chain
  --align               Enable signal-level alignment for chain evaluation
  --align-backend       Alignment backend: path-guided (default), radius, auto

Signal Processing Options (only with --graph):
  --event-pipeline      Event pipeline backend: scrappie (default), rawhash, passthrough
  --fuzzy-backend       Fuzzy quantizer backend (default: rh2)
  --fuzzy-fine-min      Fuzzy quantizer fine region min (default: -2.0)
  --fuzzy-fine-max      Fuzzy quantizer fine region max (default: 2.0)
  --fuzzy-fine-range    Fuzzy quantizer fine bin width (default: 0.4)
  --seed-k              Seed extractor k-mer size (default: 6)
  --seed-stride         Seed extractor stride (default: 1)

Debug Options:
  --dump-anchors        Dump anchors to directory (one file per read)
  --dump-chains         Dump chains to directory (one file per read)

Output Options:
  -o, --output          Output file path (format auto-detected from extension: .paf, .gaf, .gam, .json)
  --output-format       Override output format (paf, gaf, gam, json)
  --min-secondary-ratio Min chain score ratio vs primary for secondaries (default: 0.7)
```
</details>

## Simulation Workflow

Generate test reads using [squigulator](https://github.com/hasindu2008/squigulator):

```bash
# Generate simulated R9 reads from a reference
squigulator reference.fa -x dna-r9-min \
  -o reads.blow5 -n 20 -r 2000 \
  -q reads.fasta -c reads.paf \
  --ideal --seed 123

# Index and map
./piru/build/piru index -m r9.4 -o my_index reference.gfa
./piru/build/piru map -i my_index reads.blow5
```

**Key squigulator options:**
| Option | Description |
|--------|-------------|
| `-x dna-r9-min` / `-x dna-r10-min` | Chemistry profile (match with piru `--model`) |
| `-o reads.blow5` | Output signal file |
| `-n 20 -r 2000` | Generate 20 reads of ~2000 bp |
| `-q reads.fasta` | Output basecalled sequences |
| `-c reads.paf` | Output ground truth alignments |
| `--ideal` | Perfect simulation (no noise) |
| `--seed N` | Fixed seed for reproducibility |

## Documentation

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) - System architecture and module design
- [timing.md](docs/timing.md) - Profiling and timing utilities
- [threading.md](docs/threading.md) - Concurrency abstractions
- [testing.md](docs/testing.md) - Testing strategy

## Contributing

1. Run `make check` or `ctest --test-dir piru/build` before submitting changes
2. Follow `.clang-format` style (Google base, 4-space indent, 100-col limit)
3. Add tests for new functionality in `tests/`

## License

MIT
