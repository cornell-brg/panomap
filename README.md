<p align="center">
  <img src="docs/assets/logo.png" alt="PIRU" width="200">
</p>

---

<p align="center">
  <em>Real-time squiggle-to-graph classification for Nanopore adaptive sampling</em>
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

## What is PIRU?

PIRU maps raw Nanopore signal directly to pangenome graphs without basecalling.
By working at the squiggle level, it enables real-time keep/reject decisions
during adaptive sampling experiments. Supports GFA and VG graphs, SLOW5/BLOW5
signals, and R9.4/R10.4 pore chemistries.

## Quick Start

**Requirements:** CMake 3.16+, C++20 compiler, zlib, TBB (or `-DPIRU_FETCH_TBB=ON`)

```bash
# Clone and build
git clone --recursive https://github.com/xxx/piru.git
cd piru && mkdir build && cd build && cmake .. && make

# Index a graph
./piru index -m r10.4 --seed-mode path reference.gfa -o ref.pirx

# Map reads
./piru map --index ref.pirx reads.blow5 -o out.gaf

# ROI classification (adaptive sampling)
./piru annotate reference.gfa --bed targets.bed -o targets.pira
./piru map --index ref.pirx --roi targets.pira --mode enrich reads.blow5 -o out.gaf
```

**Optional dependencies:**
- VG support: `libjansson-dev`, `libhts-dev` (or `-DPIRU_ENABLE_LIBVGIO=OFF`)
- Tests: `ctest` from build directory

## Subcommands

| Command | Description |
|---------|-------------|
| `piru index` | Build .pirx index from GFA/VG graph + pore model |
| `piru map` | Map BLOW5 reads against index, optionally classify with ROI |
| `piru annotate` | Project BED target intervals onto graph -> .pira annotation |

Run `piru <command> --help` for full options.

## Simulation Workflow

Generate test reads using [squigulator](https://github.com/hasindu2008/squigulator):

```bash
squigulator reference.fa -x dna-r10-min \
  -o reads.blow5 -n 20 -r 2000 \
  -q reads.fasta -c reads.paf \
  --ideal --seed 123

./piru index -m r10.4 --seed-mode path -o ref.pirx reference.gfa
./piru map --index ref.pirx reads.blow5 -o out.gaf
```

## Documentation

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) -- System architecture and vision
- [index_format.md](docs/index_format.md) -- Binary index format specification

## Build & Test

```bash
mkdir build && cd build
cmake ..
make
ctest          # unit tests
ctest -V       # verbose
```

## Contributing

1. Run `make && ctest` in the build directory before submitting changes
2. Follow `.clang-format` style (Google base, 2-space indent, 100-col limit)
3. Add tests for new functionality in `tests/`

## License

MIT
