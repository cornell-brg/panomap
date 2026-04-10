# Integration Tests Scoreboard

Last updated: 2026-04-10

| Test | Dataset | Tokenizer | Params | Recall | Map Acc | Threshold |
|------|---------|-----------|--------|--------|---------|-----------|
| drb1_accuracy | DRB1 pangenome (1k sim) | rh2 k=8 | default | 99% | 98% | R>=90 |
| drb1_sort_chain | DRB1 pangenome (1k sim) | rh2 k=8 | --chainer sort-chain | 99% | 94% | R>=90 |
| drb1_pan_chain | DRB1 pangenome (1k sim) | rh2 k=8 | --chainer pan-chain | 99% | 93% | R>=90 |
| drb1_landmark | DRB1 pangenome (1k sim) | landmark k=4 | default | 100% | 99% | R>=90 |
| drb1_landmark_sort_chain | DRB1 pangenome (1k sim) | landmark k=4 | --chainer sort-chain | 99% | 95% | R>=90 |
| drb1_landmark_pan_chain | DRB1 pangenome (1k sim) | landmark k=4 | --chainer pan-chain | 99% | 94% | R>=90 |
| covid_accuracy | COVID linear (1k real) | rh2 k=8 | default | 94% | 88% | R>=80, A>=85 |
| covid_viral | COVID linear (1k real) | rh2 k=6 | viral chaining | 97% | 90% | R>=90, A>=85 |
| covid_landmark_viral | COVID linear (1k real) | landmark k=4 | viral + no-adaptive + diff=0.2 + no-merge + gap=1.6 | 87% | 81% | R>=85, A>=80 |

## Running

```bash
ctest -L integration         # integration tests only
ctest -L integration -V      # verbose
ctest                        # all tests (unit + integration)
```

## Notes

- DRB1 tests simulate reads with squigulator (requires squigulator in infra/tools/)
- COVID tests use pre-basecalled real r9.4 reads compared against minimap2 truth PAF
- Recall = mapped / total (how many reads did we call mapped)
- Mapping accuracy = correct / mapped (of mapped reads, how many hit the right position)
- Tests SKIP gracefully if data or tools missing
