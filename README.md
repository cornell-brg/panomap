# Panomap

<p align="left">
  <em>Nanopore signal mapping to pangenome variation graphs</em>
</p>

<p align="left">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

## What is Panomap?

Panomap maps raw Nanopore signal reads (in SLOW5/BLOW5 format) to pangenome
variation graphs without basecalling. It works at the signal level and targets
GFA graphs built with tools like PGGB and minigraph-cactus.

## Requirements

- **CMake** 3.16–3.31 (CMake 4.x works with one extra flag — see [Build options](#build-options))
- **C++20 compiler** (GCC 13+ recommended)
- **zlib**
- **oneTBB** — parallel index/map backend; uses a system install if present, otherwise auto-fetched during configure
- **libzstd** *(optional)* — only needed for zstd-compressed BLOW5 input

Developed and tested mainly on Red Hat Enterprise Linux 8 with GCC 13.3.1 and
CMake 3.26. If it doesn't build on your system, please open an issue or pull
request — we're happy to help.

## Quick Start

```bash
# Clone and build
git clone --recursive https://github.com/cornell-brg/panomap.git
cd panomap && mkdir build && cd build && cmake .. && make -j8

# Try the bundled example (SARS-CoV-2 pangenome + 20 simulated reads)
./panomap index -m r10.4 ../examples/covid/covid-pangenome.gfa -o covid.pirx
./panomap map --index covid.pirx ../examples/covid/reads.blow5 -o out.gaf

# --- with your own data ---

# Index a graph
./panomap index -m r10.4 reference.gfa -o ref.pirx

# Map reads (GAF to stdout by default)
./panomap map --index ref.pirx reads.blow5 -o out.gaf

# Map with viral preset params (good for small genomes)
./panomap map --index ref.pirx reads.blow5 \
  --chain-bw 100 --chain-max-dist 500 --chain-pen-gap 1.2 \
  --chain-pen-skip 0.3 --max-chunks 5 -o out.gaf
```

## Build options

`cmake .. && make` works out of the box. The flags below cover special cases only.

| Flag | Default | Purpose |
|------|---------|---------|
| `-DPANOMAP_USE_ZSTD=AUTO\|ON\|OFF` | `AUTO` | zstd-compressed BLOW5 support. `AUTO` enables it when libzstd is found; `ON` fails configure if libzstd is missing; `OFF` disables. |
| `-DPANOMAP_FETCH_TBB=ON\|OFF` | `ON` | Auto-fetch oneTBB when no system install is found. `OFF` requires a system oneTBB. |
| `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` | — | Required on **CMake 4.x**: the pinned oneTBB predates CMake 4 and otherwise errors with `Compatibility with CMake < 3.5 has been removed`. |

**zstd-compressed BLOW5.** BLOW5 records may be zlib- or zstd-compressed; zlib is
always supported. To read zstd-compressed files, install `libzstd-dev`
(Debian/Ubuntu) or `libzstd-devel` (RHEL/Fedora) before configuring — it is then
picked up automatically. Without it, mapping a zstd BLOW5 aborts with
`slow5lib has not been compiled with zstd support`.

## Subcommands

| Command | Description |
|---------|-------------|
| `panomap index` | Build .pirx index from GFA graph + pore model |
| `panomap map` | Map BLOW5/SLOW5 reads against index, output GAF/PAF |

Run `panomap <command> --help` for full options.

## GAF Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `pn:Z:` | string | Reference path name (or `*` if unmapped) |
| `tp:A:` | char | Alignment type: `P` primary, `S` secondary, `U` unmapped |
| `cs:i:` | int | Chain score |
| `an:i:` | int | Anchor count in chain |
| `se:f:` | float | Score per event span (chain_score / query_span) |
| `ad:f:` | float | Anchor density (anchors / ref_span) |
| `ci:f:` | float | Canonical 1D coordinate start |
| `ce:f:` | float | Canonical 1D coordinate end |
| `cc:i:` | int | 1D component ID |
| `ws:f:` | float | Weighted standout score (mapping decision confidence) |
| `nc:i:` | int | Number of competitive chains (after secondary ratio filter) |
| `ck:i:` | int | Chunks processed before decision |
| `dt:f:` | float | Processing time (seconds) |

## Simulation Workflow

Generate test reads using [squigulator](https://github.com/hasindu2008/squigulator):

```bash
squigulator reference.fa -x dna-r9-min \
  -o reads.blow5 -n 20 -r 8000 \
  --sample-rate 4000 \
  -q reads.fasta -c reads.paf \
  --seed 123

./panomap index -m r9.4 reference.gfa -o ref.pirx
./panomap map --index ref.pirx reads.blow5 -o out.gaf
```

## Build & Test

Tests use [doctest](https://github.com/doctest/doctest) (fetched automatically
during configure). Build and run them with:

```bash
mkdir build && cd build
cmake ..
make -j8
make test
```

## Reproducing paper results

The [`reproduce/`](reproduce/) directory has self-contained scripts to reproduce
the panomap results in the paper from the data on Zenodo
([10.5281/zenodo.21420009](https://doi.org/10.5281/zenodo.21420009)). See
[`reproduce/README.md`](reproduce/README.md).

## Citing Panomap

If you use Panomap in your work, please cite:

```bibtex
@article{shih2026panomap,
  title   = {Panomap: Unbiased Nanopore Signal Mapping with Pangenome Variation Graphs},
  author  = {Shih, Po Jui and Sanghani, Zephan and Guarracino, Andrea and Gamaarachchi, Hasindu and Batten, Christopher},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.07.10.737796},
  url     = {https://doi.org/10.64898/2026.07.10.737796}
}
```

## Acknowledgements

Panomap builds on code and ideas from several projects:

- [RawHash2](https://github.com/CMU-SAFARI/RawHash) -- signal event detection
  (originally from Scrappie), signal tokenization (adaptive quantization),
  chaining (originally from minimap2), and mapping result scoring.
- [minimap2](https://github.com/lh3/minimap2) -- seed index / hash table construction.
- [odgi](https://github.com/pangenome/odgi) -- PG-SGD 1D layout.

Bundled dependencies: [slow5lib](https://github.com/hasindu2008/slow5lib) (SLOW5/BLOW5
IO) and [kmer_models](https://github.com/nanoporetech/kmer_models) (ONT pore models).
Example and test reads are simulated with
[squigulator](https://github.com/hasindu2008/squigulator).

## License

MIT
