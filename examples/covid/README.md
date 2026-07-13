# COVID pangenome example

A small, self-contained example you can run right after building.

| File | Description |
|------|-------------|
| `covid-pangenome.gfa` | SARS-CoV-2 variation graph over 6 WHO lineages (Wuhan, Alpha, Delta, BA.1, BA.5, XBB.1.5), built with PGGB. 566 nodes. |
| `reads.blow5` | 20 simulated R10 signal reads (2 kb each) drawn from the 6 lineages. |

## Run it

```bash
# from the build directory
./panomap index -m r10.4 ../examples/covid/covid-pangenome.gfa -o covid.pirx
./panomap map --index covid.pirx ../examples/covid/reads.blow5 -o out.gaf
```

All 20 reads map to their source lineage path in the graph.

## Regenerating the reads

The reads were simulated with [squigulator](https://github.com/hasindu2008/squigulator)
from the 6 member genomes:

```bash
squigulator members.fa -x dna-r10-min \
  -o reads.blow5 -n 20 -r 2000 --ideal --seed 123
```
