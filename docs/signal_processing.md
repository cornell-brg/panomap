# Signal Processing

Notes on the signal-processing components shared by indexing and mapping.

## Seed Extractor

```
FuzzyQuantizedSignal     EventSeries*
 (int16_t tokens)         (optional)
        │                     │
        └──────────┬──────────┘
                   ▼
            ┌──────────────┐
            │     Seed     │
            │  Extractor   │
            │   (k-mer)    │
            └──────────────┘
                config:
              k, stride,
              qbits, window
                   │
                   ▼
              SeedBuffer
           (vector of Seeds)
         {hash, position, span}
              per window
```

- **Purpose**: convert fuzzy-quantized signal tokens into hashable seeds for lookup. Runs after fuzzy quantization and before SeedStore interaction.
- **Inputs**:
  - `FuzzyQuantizedSignal`: `std::vector<int16_t> tokens` produced by the fuzzy quantizer.
  - `EventSeries*` (optional): event boundaries (`start`, `length`, `mean`, `stdv`). Current k-mer backend ignores it; future event-aware extractors may use it.
- **Outputs**:
  - `SeedBuffer`: `std::vector<Seed>` with each `Seed` carrying `{hash, position, span}`.
- **Interface**: `SeedExtractor` (`include/signal/seed_extractors/seed_extractor.hpp`) with `extract(...)`, `config()`, and `name()`. Config fields include backend name, `k`, `stride`, `qbits`, optional window/minimizer params.
- **Default backend (k-mer)**:
  - Slides a window of length `k` over `tokens` with step `stride` (>=1).
  - **Dual packing strategy** (to avoid UB from shift-by-64+):
    - **Normal case** (k*qbits < 64): bit-packs tokens via `(packed << qbits) | token`.
    - **Overflow case** (k*qbits >= 64 or qbits==0): incrementally mixes tokens via `hash64(packed ^ token, mask)`.
  - Final hash uses legacy `hash64` mix masked to 32 bits (stored in 64-bit `Seed.hash` field for future backend compatibility).
  - Emits one seed per window: `hash`, `position` (window start), `span` (`k`).
  - Returns empty if `k==0` or not enough tokens.
- **Factory**: `make_seed_extractor(config)` (in `signal/seed_extractors/seed_extractor_factory.cpp`) currently always constructs the k-mer backend. Unknown backends log a warning and fall back to k-mer.
