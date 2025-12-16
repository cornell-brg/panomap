#!/usr/bin/env python3
"""
Signal alignment visualization tool for debugging PIRU mapping.

Plots node signals from GFA dumps vs read signals from signal_dump_*.txt files.

Usage:
    # Compare node vs read signals
    python plot_signal_alignment.py --node-id 31 \
           --graph-dir ../build/graph_dumps \
           --signal-dump ../build/signal_dump_0.txt \
           --read-repr normalized \
           --node-repr fuzzy

    # Plot just read signal progression
    python plot_signal_alignment.py --signal-dump ../build/signal_dump_0.txt

    # Save to file
    python plot_signal_alignment.py --node-id 31 \
           --graph-dir ../build/graph_dumps \
           --signal-dump ../build/signal_dump_0.txt \
           --read-repr raw \
           --node-repr raw \
           --output comparison.png
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class SignalDumpLoader:
    """Load signals from signal_dump_*.txt files."""

    def __init__(self, dump_path: Path):
        self.dump_path = dump_path
        self.read_id = None
        self.raw_signal = None
        self.event_means = None
        self.normalized = None
        self.fuzzy_tokens = None
        self._load()

    def _load(self):
        """Load all signal representations from dump file."""
        with open(self.dump_path) as f:
            lines = f.readlines()

        if len(lines) < 5:
            raise ValueError(f"Signal dump file has only {len(lines)} lines, expected 5")

        # Line 1: Read ID
        self.read_id = lines[0].strip()

        # Line 2: Raw signal (picoamps)
        self.raw_signal = np.array([float(x) for x in lines[1].strip().split(',')])

        # Line 3: Event means (picoamps)
        self.event_means = np.array([float(x) for x in lines[2].strip().split(',')])

        # Line 4: Normalized signal (z-scores)
        self.normalized = np.array([float(x) for x in lines[3].strip().split(',')])

        # Line 5: Fuzzy quantized tokens (integer IDs)
        self.fuzzy_tokens = np.array([int(x) for x in lines[4].strip().split(',')])

    def get_signal(self, repr_type: str) -> np.ndarray:
        """Get signal by representation type."""
        if repr_type == 'raw':
            return self.raw_signal
        elif repr_type == 'events':
            return self.event_means
        elif repr_type == 'normalized':
            return self.normalized
        elif repr_type == 'fuzzy':
            return self.fuzzy_tokens.astype(float)
        else:
            raise ValueError(f"Unknown representation: {repr_type}")

    def get_label(self, repr_type: str) -> str:
        """Get axis label for representation type."""
        labels = {
            'raw': 'Signal (pA)',
            'events': 'Event Mean (pA)',
            'normalized': 'Normalized Signal (z-score)',
            'fuzzy': 'Fuzzy Token ID'
        }
        return labels.get(repr_type, repr_type)


class GFASignalLoader:
    """Load signal annotations from GFA dumps."""

    @staticmethod
    def load_node_sequence(gfa_path: Path, node_id: int) -> Optional[str]:
        """Extract node sequence from GFA."""
        with open(gfa_path) as f:
            for line in f:
                if line.startswith('S\t'):
                    parts = line.strip().split('\t')
                    if int(parts[1]) == node_id:
                        return parts[2]
        return None

    @staticmethod
    def load_node_signal(gfa_path: Path, node_id: int, signal_tag: str = 'SG') -> Optional[np.ndarray]:
        """
        Extract signal annotation from GFA node.

        Signal tags:
        - SG:Z: Raw signal (picoamps, comma-separated)
        - FZ:Z: Fuzzy quantized tokens (integer IDs, comma-separated)
        - AQ:Z: Alignment quantized signal (picoamps, comma-separated)
        """
        with open(gfa_path) as f:
            for line in f:
                if line.startswith('S\t'):
                    parts = line.strip().split('\t')
                    if int(parts[1]) == node_id:
                        # Look for signal tag
                        for tag in parts[3:]:
                            if tag.startswith(f'{signal_tag}:'):
                                signal_str = tag.split(':', 2)[2]
                                return np.array([float(x) for x in signal_str.split(',')])
        return None


class SignalPlotter:
    """Visualize node and read signals for debugging."""

    def __init__(self, graph_dir: Optional[Path] = None, signal_dump: Optional[Path] = None):
        self.graph_dir = graph_dir
        self.signal_dump = signal_dump

        # GFA dump paths (from PIRU_DUMP_GRAPHS=ON)
        if graph_dir:
            self.imported_gfa = graph_dir / "imported_graph.gfa"
            self.raw_signals_gfa = graph_dir / "raw_signals.gfa"
            self.fuzzy_gfa = graph_dir / "fuzzy_quantized.gfa"
            self.aln_gfa = graph_dir / "aln_quantized.gfa"

        # Load signal dump if provided
        self.dump_loader = None
        if signal_dump:
            self.dump_loader = SignalDumpLoader(signal_dump)

    def plot_read_progression(self, figsize=(14, 10)):
        """Plot all read signal processing stages."""
        if not self.dump_loader:
            print("Error: No signal dump file provided", file=sys.stderr)
            return None

        fig = plt.figure(figsize=figsize)
        gs = GridSpec(4, 1, figure=fig, hspace=0.3)

        # 1. Raw signal
        ax1 = fig.add_subplot(gs[0])
        ax1.plot(self.dump_loader.raw_signal, linewidth=0.5, color='black')
        ax1.set_title(f"Read {self.dump_loader.read_id} - Raw Signal (n={len(self.dump_loader.raw_signal)})")
        ax1.set_ylabel("Signal (pA)")
        ax1.grid(True, alpha=0.3)

        # 2. Event means
        ax2 = fig.add_subplot(gs[1])
        ax2.plot(self.dump_loader.event_means, linewidth=1, color='blue', marker='.')
        ax2.set_title(f"Event Means (n={len(self.dump_loader.event_means)})")
        ax2.set_ylabel("Event Mean (pA)")
        ax2.grid(True, alpha=0.3)

        # 3. Normalized signal
        ax3 = fig.add_subplot(gs[2])
        ax3.plot(self.dump_loader.normalized, linewidth=0.5, color='green')
        ax3.set_title(f"Normalized Signal (n={len(self.dump_loader.normalized)})")
        ax3.set_ylabel("Normalized (z-score)")
        ax3.grid(True, alpha=0.3)

        # 4. Fuzzy tokens
        ax4 = fig.add_subplot(gs[3])
        ax4.step(range(len(self.dump_loader.fuzzy_tokens)), self.dump_loader.fuzzy_tokens,
                where='mid', color='orange')
        ax4.set_title(f"Fuzzy Quantized Tokens (n={len(self.dump_loader.fuzzy_tokens)})")
        ax4.set_ylabel("Token ID")
        ax4.set_xlabel("Sample Index")
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def plot_node_vs_read(self, node_id: int, read_repr: str = 'raw', node_repr: str = 'raw',
                         window_start: Optional[int] = None, window_end: Optional[int] = None,
                         figsize=(14, 8)):
        """
        Compare node signal vs read signal.

        Args:
            node_id: Graph node ID
            read_repr: Read representation ('raw', 'events', 'normalized', 'fuzzy')
            node_repr: Node representation ('raw', 'fuzzy', 'alignment')
            window_start: Start sample index (None = full signal)
            window_end: End sample index (None = full signal)
        """
        if not self.graph_dir:
            print("Error: No graph directory provided", file=sys.stderr)
            return None
        if not self.dump_loader:
            print("Error: No signal dump file provided", file=sys.stderr)
            return None

        # Load node signal based on representation
        node_signal = None
        node_seq = None
        if node_repr == 'raw':
            node_signal = GFASignalLoader.load_node_signal(self.raw_signals_gfa, node_id, 'SG')
            node_seq = GFASignalLoader.load_node_sequence(self.imported_gfa, node_id)
        elif node_repr == 'fuzzy':
            node_signal = GFASignalLoader.load_node_signal(self.fuzzy_gfa, node_id, 'FZ')
            node_seq = GFASignalLoader.load_node_sequence(self.imported_gfa, node_id)
        elif node_repr == 'alignment':
            node_signal = GFASignalLoader.load_node_signal(self.aln_gfa, node_id, 'AQ')
            node_seq = GFASignalLoader.load_node_sequence(self.imported_gfa, node_id)

        if node_signal is None:
            print(f"Error: Could not load node {node_id} signal (repr={node_repr})", file=sys.stderr)
            return None

        # Load read signal
        read_signal = self.dump_loader.get_signal(read_repr)

        # Apply window if specified
        if window_start is not None or window_end is not None:
            start = window_start if window_start is not None else 0
            end = window_end if window_end is not None else len(read_signal)
            read_signal = read_signal[start:end]
            read_x = np.arange(start, end)
        else:
            read_x = np.arange(len(read_signal))

        # Create figure
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize, height_ratios=[1, 1, 1.5])

        # Top: Node signal
        node_x = np.arange(len(node_signal))
        ax1.plot(node_x, node_signal, linewidth=1.5, color='red')
        title = f"Node {node_id} - {node_repr} signal"
        if node_seq:
            title += f" (seq: {node_seq[:30]}{'...' if len(node_seq) > 30 else ''})"
        ax1.set_title(title)
        ax1.set_ylabel(self.dump_loader.get_label(node_repr) if node_repr == read_repr else node_repr)
        ax1.grid(True, alpha=0.3)

        # Middle: Read signal
        ax2.plot(read_x, read_signal, linewidth=1.5, color='blue')
        ax2.set_title(f"Read {self.dump_loader.read_id} - {read_repr} signal")
        ax2.set_ylabel(self.dump_loader.get_label(read_repr))
        ax2.grid(True, alpha=0.3)

        # Bottom: Normalized overlay
        # Normalize both to [0, 1] for visual comparison
        node_norm = (node_signal - np.min(node_signal)) / (np.max(node_signal) - np.min(node_signal) + 1e-6)
        read_norm = (read_signal - np.min(read_signal)) / (np.max(read_signal) - np.min(read_signal) + 1e-6)

        ax3.plot(node_x, node_norm, linewidth=1.5, color='red', alpha=0.7, label=f'Node {node_id} (normalized)')

        # Align read to node length if different
        if len(read_signal) == len(node_signal):
            ax3.plot(read_x, read_norm, linewidth=1.5, color='blue', alpha=0.7, label='Read (normalized)')
            # Compute correlation
            corr = np.corrcoef(node_norm, read_norm)[0, 1]
            ax3.text(0.02, 0.98, f"Correlation: {corr:.3f}",
                    transform=ax3.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        else:
            # Different lengths - plot separately
            read_x_scaled = np.linspace(0, len(node_signal)-1, len(read_signal))
            ax3.plot(read_x_scaled, read_norm, linewidth=1.5, color='blue', alpha=0.7, label='Read (normalized, scaled)')
            ax3.text(0.02, 0.98, f"Note: Different lengths (node={len(node_signal)}, read={len(read_signal)})",
                    transform=ax3.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

        ax3.set_title(f"Normalized Overlay Comparison")
        ax3.set_ylabel("Normalized Signal [0,1]")
        ax3.set_xlabel("Sample Index")
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig


def main():
    parser = argparse.ArgumentParser(
        description="Visualize PIRU signal processing for debugging",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot read signal progression
  python plot_signal_alignment.py --signal-dump signal_dump_0.txt

  # Compare node raw signal vs read normalized signal
  python plot_signal_alignment.py --node-id 31 \\
      --graph-dir build/graph_dumps \\
      --signal-dump signal_dump_0.txt \\
      --read-repr normalized --node-repr raw

  # Compare fuzzy quantized signals
  python plot_signal_alignment.py --node-id 31 \\
      --graph-dir build/graph_dumps \\
      --signal-dump signal_dump_0.txt \\
      --read-repr fuzzy --node-repr fuzzy

  # Focus on specific window of read signal
  python plot_signal_alignment.py --node-id 31 \\
      --graph-dir build/graph_dumps \\
      --signal-dump signal_dump_0.txt \\
      --read-repr raw --node-repr raw \\
      --window-start 0 --window-end 500
        """
    )

    parser.add_argument('--signal-dump', type=Path, required=True,
                       help="Path to signal_dump_*.txt file")
    parser.add_argument('--node-id', type=int, help="Graph node ID to compare")
    parser.add_argument('--graph-dir', type=Path,
                       help="Directory containing GFA dumps (required with --node-id)")
    parser.add_argument('--read-repr', choices=['raw', 'events', 'normalized', 'fuzzy'],
                       default='raw', help="Read signal representation (default: raw)")
    parser.add_argument('--node-repr', choices=['raw', 'fuzzy', 'alignment'],
                       default='raw', help="Node signal representation (default: raw)")
    parser.add_argument('--window-start', type=int, help="Start sample index for read window")
    parser.add_argument('--window-end', type=int, help="End sample index for read window")
    parser.add_argument('--output', type=Path, help="Save plot to file (default: show interactive)")

    args = parser.parse_args()

    # Validate arguments
    if args.node_id is not None and args.graph_dir is None:
        parser.error("--graph-dir is required when --node-id is specified")

    plotter = SignalPlotter(args.graph_dir, args.signal_dump)

    # Plot based on what's requested
    if args.node_id is not None:
        # Node vs read comparison
        print(f"Plotting node {args.node_id} ({args.node_repr}) vs read ({args.read_repr})...")
        fig = plotter.plot_node_vs_read(
            args.node_id,
            read_repr=args.read_repr,
            node_repr=args.node_repr,
            window_start=args.window_start,
            window_end=args.window_end
        )
    else:
        # Just read progression
        print(f"Plotting read signal progression...")
        fig = plotter.plot_read_progression()

    if fig is None:
        return 1

    # Save or show
    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {args.output}")
    else:
        plt.show()

    return 0


if __name__ == '__main__':
    sys.exit(main())
