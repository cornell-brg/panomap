# Test Reads for PIRU Quick Testing

This directory contains several small sets of simulated nanopore reads that make
it easy to sanity-check PIRU’s mapper without running full experiments. Each set
consists of matching `.blow5`, `.fasta`, and `.paf` files produced by
`squigulator` from the same reference (`../DRB1-3123.fa`).

## Available datasets

| Prefix | Files | Chemistry profile | Target read length | Reads | Notes |
| --- | --- | --- | --- | --- | --- |
| `quick_r9_2k` | `quick_r9_2k.{blow5,fasta,paf}` | `dna-r9-min` | ~2000 bp | 5 | Historic default used in CI/smoke tests |
| `quick_r9_1k` | `quick_r9_1k.{blow5,fasta,paf}` | `dna-r9-min` | ~1000 bp | 5 | Shorter reads for stress-testing normalization |
| `quick_r9_4k` | `quick_r9_4k.{blow5,fasta,paf}` | `dna-r9-min` | ~4000 bp | 5 | Longer R9 reads |
| `quick_r10_1k` | `quick_r10_1k.{blow5,fasta,paf}` | `dna-r10-min` | ~1000 bp | 5 | Short R10 chemistry reads |
| `quick_r10_4k` | `quick_r10_4k.{blow5,fasta,paf}` | `dna-r10-min` | ~4000 bp | 5 | Longer R10 chemistry reads |

Pick any dataset whose chemistry matches the `--pore-model` flag you plan to
use (`r9` or `r10`). All sets were generated in **ideal** mode so their expected
mapping accuracy is 100%.

## Quick test examples

From the PIRU build directory:

```bash
module load slow5lib           # ensure libslow5 is available

# Example: run against the default R9 dataset
./piru map --pore-model r9 \
  ../data/HLA/DRB1-3123.gfa \
  ../data/HLA/test_reads/quick_r9_2k.blow5

# Example: run against the R10 4kb dataset
./piru map --pore-model r10 \
  ../data/HLA/DRB1-3123.gfa \
  ../data/HLA/test_reads/quick_r10_4k.blow5
```

All 5 reads should map successfully in under a second on a typical workstation.

## Generation log

All datasets were generated with **squigulator 0.4.0** on 2025-11-20 using the
commands below (run from `data/HLA/test_reads/` with `module load
squigulator/0.4.0`). Seeds are fixed for reproducibility.

```bash
# R9 chemistries
squigulator ../DRB1-3123.fa -x dna-r9-min \
  -o quick_r9_2k.blow5 -n 5 -r 2000 \
  -q quick_r9_2k.fasta -c quick_r9_2k.paf \
  --ideal --seed 41

squigulator ../DRB1-3123.fa -x dna-r9-min \
  -o quick_r9_1k.blow5 -n 5 -r 1000 \
  -q quick_r9_1k.fasta -c quick_r9_1k.paf \
  --ideal --seed 42

squigulator ../DRB1-3123.fa -x dna-r9-min \
  -o quick_r9_4k.blow5 -n 5 -r 4000 \
  -q quick_r9_4k.fasta -c quick_r9_4k.paf \
  --ideal --seed 43

# R10 chemistries
squigulator ../DRB1-3123.fa -x dna-r10-min \
  -o quick_r10_1k.blow5 -n 5 -r 1000 \
  -q quick_r10_1k.fasta -c quick_r10_1k.paf \
  --ideal --seed 44

squigulator ../DRB1-3123.fa -x dna-r10-min \
  -o quick_r10_4k.blow5 -n 5 -r 4000 \
  -q quick_r10_4k.fasta -c quick_r10_4k.paf \
  --ideal --seed 45
```

Re-run the relevant command to refresh any dataset (adjust `-n`, `-r`, or seeds
as needed). Remember to match the `--pore-model` flag to the chemistry profile
when mapping these reads.
