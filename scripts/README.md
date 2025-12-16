# PIRU Scripts

Utility scripts for debugging and analysis.

## plot_signal_alignment.py

Visualization tool for debugging signal processing and mapping.

### Purpose

Helps diagnose signal processing issues by plotting:
- Node signals at different pipeline stages (raw → fuzzy quantized → alignment quantized)
- Read signals from BLOW5/SLOW5 files
- Anchor matches overlaid to verify signal similarity

### Prerequisites

```bash
pip install matplotlib numpy pyslow5
```

### Enable Graph Dumping

PIRU must be built with `PIRU_DUMP_GRAPHS=ON` to generate GFA annotations:

```bash
cd build
cmake -DPIRU_DUMP_GRAPHS=ON ..
cmake --build .
```

When enabled, indexing and mapping will output GFA files to `build/graph_dumps/`:
- `imported_graph.gfa` - Original graph structure
- `transformed_graph.gfa` - After VG transform (if applicable)
- `raw_signals.gfa` - Node signals after squigglization
- `fuzzy_quantized.gfa` - Fuzzy quantized tokens (seed extraction)
- `aln_quantized.gfa` - Alignment quantized signals

### Usage

**1. Plot node signal processing stages:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --graph-dir build/graph_dumps
```

Shows:
- Node sequence
- Raw signal (picoamps)
- Fuzzy quantized tokens
- Alignment quantized signal

**2. Plot read signal:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --read-id "S1_1!gi|568815592:32578768-32589835!54!1312!-" \
    --reads tests/data/HLA/test_reads/quick_r9_1k.blow5 \
    --graph-dir build/graph_dumps
```

Shows:
- Read raw ADC signal
- Converted to picoamps

**3. Compare node vs read at anchor position:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --read-id "S1_1!gi|568815592:32578768-32589835!54!1312!-" \
    --reads tests/data/HLA/test_reads/quick_r9_1k.blow5 \
    --graph-dir build/graph_dumps \
    --anchor-query 482 \
    --anchor-ref 1173 \
    --window-size 100 \
    --output anchor_comparison.png
```

Shows:
- Full node and read signals
- Overlaid normalized signals at anchor window
- Correlation metric
- Fuzzy quantized tokens

### Debugging Workflow

1. **Identify suspicious anchor** from mapping output:
   ```
   S1_1!...  anchors_detail=[q=482,r=1173,l=6,p=0]
   ```

2. **Find corresponding node ID** for ref_coord=1173 on path=0:
   - Check linearization coordinates (future: add query tool)
   - Or inspect GFA dumps directly

3. **Visualize anchor match**:
   ```bash
   python scripts/plot_signal_alignment.py \
       --node-id <node_id> \
       --read-id "S1_1!..." \
       --reads <reads.blow5> \
       --graph-dir build/graph_dumps \
       --anchor-query 482 \
       --anchor-ref 1173
   ```

4. **Check for issues**:
   - Low correlation (<0.5): anchor match is spurious
   - Identical fuzzy tokens across all nodes: quantization collapse
   - Flat signals: event detection or normalization broken
   - Mismatched signal scales: calibration issue

### Common Issues

**"pyslow5 not available"**
```bash
pip install pyslow5
```

**"Could not load read"**
- Check read ID matches exactly (case-sensitive, special chars)
- Verify BLOW5 file path is correct

**"No GFA dumps found"**
- Rebuild with `cmake -DPIRU_DUMP_GRAPHS=ON`
- Run index or map command to generate dumps
- Check `build/graph_dumps/` directory exists

### Output Formats

- Interactive: omit `--output` to show plots with matplotlib GUI
- PNG: `--output plot.png`
- Use `.pdf` or `.svg` for vector output
