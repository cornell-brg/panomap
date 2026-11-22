# Panomap

<p align="left">
  <em>Nanopore signal mapping to pangenome variation graphs</em>
</p>

<p align="left">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

Panomap maps raw Nanopore signal reads (SLOW5/BLOW5) to pangenome variation graphs
(GFA) without basecalling. It works on graphs built with tools such as PGGB and
minigraph-cactus.

*Panomap is proof-of-concept research software. It does not support human-scale
pangenomes yet.*

## Install

Each command in this README starts with `%`, the shell prompt. To paste commands into
a shell, run `alias %=""` first.

You need:

- CMake 3.16 to 3.31 (CMake 4.x needs one extra flag, see below)
- A C++20 compiler (GCC 13 or newer recommended)
- zlib
- oneTBB (configure fetches it if there is no system install)
- libzstd, only to read zstd-compressed BLOW5

Build and test:

```bash
% git clone --recursive https://github.com/cornell-brg/panomap.git
% cd panomap && mkdir build && cd build
% cmake .. && make -j8
% make test
```

The tests use [doctest](https://github.com/doctest/doctest), which configure fetches.
Panomap is developed on Red Hat Enterprise Linux 8 with GCC 13.3.1 and CMake 3.26. If it
does not build on your system, open an issue.

Optional CMake flags:

| Flag | Default | Purpose |
|------|---------|---------|
| `-DPANOMAP_USE_ZSTD=AUTO\|ON\|OFF` | `AUTO` | zstd-compressed BLOW5. With `AUTO`, a found libzstd turns it on. Without libzstd, `ON` fails configure. |
| `-DPANOMAP_FETCH_TBB=ON\|OFF` | `ON` | Without a system oneTBB, fetch it. |
| `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` | none | Needed on CMake 4.x, because the pinned oneTBB predates CMake 4. |

Without libzstd, mapping a zstd-compressed BLOW5 stops with `slow5lib has not been
compiled with zstd support`. Install `libzstd-dev` (Debian, Ubuntu) or `libzstd-devel`
(RHEL, Fedora) and configure again.

## Quick start

Run the bundled example, a SARS-CoV-2 pangenome with 20 simulated reads:

```bash
% ./panomap index -m r10.4 ../examples/covid/covid-pangenome.gfa -o covid.pirx
% ./panomap map --index covid.pirx ../examples/covid/reads.blow5 -o out.gaf
```

For your own data, build the index with a preset. `map` then uses the same preset.

```bash
% ./panomap index -x r10 graph.gfa -o graph.pirx
% ./panomap map --index graph.pirx reads.blow5 -o out.paf
```

## Usage

### Commands

| Command | What it does |
|---------|--------------|
| `panomap index` | Builds a `.pirx` index from a GFA graph and a pore model |
| `panomap map` | Maps SLOW5/BLOW5 reads against an index and writes PAF or GAF |
| `panomap inspect` | Shows an index's settings and section sizes. With `--compare`, checks two indexes section by section. |

Run `panomap <command> --help` for all options.

### Presets

A preset sets the pore model, the seeds, and the map settings together. `index -x
<preset>` stores the preset in the index, and `map` uses it. A preset named at `map`
time replaces it, and flags after a preset override it.

| Preset | Chemistry | Use |
|--------|-----------|-----|
| `r9` | R9.4.1 | Genomes of bacterial size and up |
| `r9-fast` | R9.4.1 | About 3x faster mapping. Best with 4 or more chunks per read. |
| `r9-small` | R9.4.1 | Small genomes, such as viruses |
| `r10` | R10.4.1 | Genomes of bacterial size |
| `r10-fast` | R10.4.1 | About 5x faster mapping, at lower recall |
| `r10-large` | R10.4.1 | Large genomes |

### Output

For an output name that ends in `.paf`, `map` writes PAF. For any other name, it writes
GAF. Without `-o`, it writes GAF to stdout.

- Both formats give the read length in signal samples (column 2) and the read start
  and end as event positions (columns 3 and 4).
- PAF gives the chosen path as the target, with positions on that path.
- GAF gives the walk of every path step that the chain covers (column 6) and the
  positions on that walk (columns 7 to 9). The `pn:Z` tag names the path.
- `map --gaf-anchor-nodes` makes GAF list only the nodes that hold anchors, with path
  positions in columns 7 to 9, as in PAF.
- Unmapped reads get a line with `tp:A:U`.

| Tag | Type | Meaning |
|-----|------|---------|
| `tp` | A | `P` primary, `S` secondary (GAF only), `U` unmapped |
| `pn` | Z | Chosen path (GAF only). `*` for an unmapped read. |
| `s1` | i | Chain score |
| `an` | i | Anchors in the chain |
| `se` | f | Chain score per read event |
| `ad` | f | Anchors per reference base |
| `ws` | f | Weighted standout of the best chain |
| `nc` | i | Chains kept for the read |
| `df` | A | Decision: `G` single-chain gate, `S` standout, `F` fallback, `U` unmapped |
| `mt` | f | Map time of the read, in ms |
| `ci` | i | Chunks processed |
| `ls`, `le`, `lc` | f, f, i | 1D layout start, end, and component (GAF, for an index with a layout) |

The tags `ws` to `ci` appear once per read, on its primary or unmapped line. `mt` and
`ci` use the names that RawHash2 and its scripts read.

### Repeat builds

Index builds with more than one thread differ slightly in the 1D layout. The layout
uses path-guided SGD with lock-free parallel updates, as in odgi. Mapping uses the
layout to merge duplicate chains, so a few reads can map differently between two
builds. Builds with one thread (`-t 1`) are identical. `panomap inspect --compare
a.pirx b.pirx` shows which parts of two indexes differ.

## Simulating reads

Simulate reads with [squigulator](https://github.com/hasindu2008/squigulator), then map
them:

```bash
% squigulator reference.fa -x dna-r9-min -o reads.blow5 -n 20 -r 8000 \
    --sample-rate 4000 -q reads.fasta -c reads.paf --seed 123

% ./panomap index -x r9 reference.gfa -o ref.pirx
% ./panomap map --index ref.pirx reads.blow5 -o out.paf
```

## Reproducing paper results

The [`reproduce/`](reproduce/) directory has scripts that reproduce the Panomap results
in the paper from the data on Zenodo
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

- [RawHash2](https://github.com/CMU-SAFARI/RawHash): signal event detection (originally
  from Scrappie), signal tokenization (adaptive quantization), chaining (originally
  from minimap2), and mapping result scoring.
- [minimap2](https://github.com/lh3/minimap2): seed index and hash table construction.
- [odgi](https://github.com/pangenome/odgi): PG-SGD 1D layout.

Bundled dependencies: [slow5lib](https://github.com/hasindu2008/slow5lib) (SLOW5/BLOW5
I/O) and [kmer_models](https://github.com/nanoporetech/kmer_models) (ONT pore models).
Example and test reads are simulated with
[squigulator](https://github.com/hasindu2008/squigulator).

## License

MIT
