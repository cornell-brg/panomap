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

PIRU maps raw Nanopore signal directly to pangenome graphs without basecalling.
By working at the squiggle level, it enables real-time keep/reject decisions
during adaptive sampling experiments. Supports GFA graphs, SLOW5/BLOW5 signals,
and R9.4/R10.4 pore chemistries.

## Quick Start

**Requirements:** CMake 3.16+, C++20 compiler, zlib, TBB (or `-DPIRU_FETCH_TBB=ON`)

```bash
# Clone and build
git clone --recursive https://github.com/xxx/piru.git
cd piru && mkdir build && cd build && cmake .. && make -j8

# Index a graph (RH2 tokenizer, default)
./piru index -m r10.4 reference.gfa -o ref.pirx

# Index with landmark tokenizer (amplitude-based peak seeding)
./piru index -m r9.4 --tokenizer landmark --seed-k 4 reference.gfa -o ref.pirx

# Map reads (GAF to stdout by default)
./piru map --index ref.pirx reads.blow5 -o out.gaf
```

## Subcommands

| Command | Description |
|---------|-------------|
| `piru index` | Build .pirx index from GFA graph + pore model |
| `piru map` | Map BLOW5/SLOW5 reads against index, output GAF/PAF |

Run `piru <command> --help` for full options.

## Tokenizers

PIRU supports two tokenization strategies for converting signal to seeds:

| Tokenizer | Description | Seed bits | Typical use |
|-----------|-------------|-----------|-------------|
| `rh2` (default) | Adaptive quantization of event values | 24 (k=6, 4-bit tokens) | General purpose |
| `landmark` | Log-quantized peak rise/fall amplitudes | 16 (k=4, 4-bit tokens) | Fewer seeds, faster chaining |

The landmark tokenizer detects valleys in the signal, extracts peaks between
them, and encodes each peak by its rise and fall amplitudes. This produces
~3x fewer index seeds and ~4x fewer hits per read compared to RH2, with
equivalent or better mapping accuracy.

## Simulation Workflow

Generate test reads using [squigulator](https://github.com/hasindu2008/squigulator):

```bash
squigulator reference.fa -x dna-r10-min \
  -o reads.blow5 -n 20 -r 2000 \
  -q reads.fasta -c reads.paf \
  --ideal --seed 123

./piru index -m r10.4 reference.gfa -o ref.pirx
./piru map --index ref.pirx reads.blow5 -o out.gaf
```

## Documentation

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) -- System architecture
- [index_format.md](docs/index_format.md) -- Binary index format specification

## Build & Test

```bash
mkdir build && cd build
cmake ..
make -j8
ctest          # unit tests
ctest -V       # verbose
```

## Contributing

1. Run `make && ctest` in the build directory before submitting changes
2. Follow `.clang-format` style (Google base, 2-space indent, 100-col limit)
3. Add tests for new functionality in `tests/`

## License

MIT
