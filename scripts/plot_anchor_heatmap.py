#!/usr/bin/env python3
"""Plot heatmap of anchor positions across haplotype paths.

This tool visualizes where read anchors map across different haplotype paths,
helping diagnose whether anchors are properly distributed or concentrated on
a single path.

Supports grouped visualization of forward (+) and reverse (-) strands per haplotype.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# Set font to Fira Sans Condensed
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Fira Sans Condensed', 'DejaVu Sans']


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class Path:
    """Represents a haplotype path in the graph."""
    path_id: int
    name: str
    length: int  # Total linear length in base pairs
    strand: str = "+"  # "+" for forward, "-" for reverse
    base_name: str = ""  # Name without _reverse suffix

    def __post_init__(self):
        if self.name.endswith("_reverse"):
            self.strand = "-"
            self.base_name = self.name[:-8]  # Remove "_reverse"
        else:
            self.strand = "+"
            self.base_name = self.name


@dataclass
class Anchor:
    """Represents a single anchor mapping between read and reference."""
    read_id: str
    path_id: int
    path_name: str
    query_pos: int      # Position in read
    ref_coord: int      # Linear coordinate on reference path
    length: int         # Anchor span
    node_id: int = 0    # Graph node ID (optional)


@dataclass
class HaplotypeGroup:
    """Groups forward and reverse paths for a single haplotype."""
    base_name: str
    forward_path: Optional[Path] = None
    reverse_path: Optional[Path] = None

    @property
    def length(self) -> int:
        """Return the length (should be same for both strands)."""
        if self.forward_path:
            return self.forward_path.length
        if self.reverse_path:
            return self.reverse_path.length
        return 0


# =============================================================================
# Input Parsing
# =============================================================================

def parse_anchor_dump(filepath: str) -> Tuple[List[Path], List[Anchor]]:
    """Parse anchor dump TSV file from piru.

    Format:
        #PATH	path_id	path_name	length
        #ANCHORS	read_id	path_id	path_name	query_pos	ref_coord	length	node_id
        <anchor data rows>

    Returns:
        Tuple of (paths, anchors)
    """
    paths = []
    anchors = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#PATH\t"):
                parts = line.split("\t")
                # #PATH	path_id	path_name	length
                path_id = int(parts[1])
                path_name = parts[2]
                length = int(parts[3])
                paths.append(Path(path_id=path_id, name=path_name, length=length))

            elif line.startswith("#ANCHORS"):
                # Header line, skip
                continue

            elif not line.startswith("#"):
                # Anchor data row
                parts = line.split("\t")
                if len(parts) >= 7:
                    anchors.append(Anchor(
                        read_id=parts[0],
                        path_id=int(parts[1]),
                        path_name=parts[2],
                        query_pos=int(parts[3]),
                        ref_coord=int(parts[4]),
                        length=int(parts[5]),
                        node_id=int(parts[6]) if len(parts) > 6 else 0
                    ))

    return paths, anchors


def group_paths_by_haplotype(paths: List[Path]) -> List[HaplotypeGroup]:
    """Group forward and reverse paths by base haplotype name."""
    groups: Dict[str, HaplotypeGroup] = {}

    for path in paths:
        if path.base_name not in groups:
            groups[path.base_name] = HaplotypeGroup(base_name=path.base_name)

        if path.strand == "+":
            groups[path.base_name].forward_path = path
        else:
            groups[path.base_name].reverse_path = path

    # Sort by base name for consistent ordering
    return sorted(groups.values(), key=lambda g: g.base_name)


# =============================================================================
# Test Data Generation
# =============================================================================

def generate_test_paths() -> List[Path]:
    """Generate test haplotype paths with forward and reverse."""
    paths = []
    base_names = [
        ("HLA-DRB1*01:01", 15000),
        ("HLA-DRB1*03:01", 14800),
        ("HLA-DRB1*04:01", 15200),
        ("HLA-DRB1*07:01", 14900),
        ("HLA-DRB1*15:01", 15100),
    ]

    path_id = 0
    for name, length in base_names:
        # Forward
        paths.append(Path(path_id=path_id, name=name, length=length))
        path_id += 1

    for name, length in base_names:
        # Reverse
        paths.append(Path(path_id=path_id, name=f"{name}_reverse", length=length))
        path_id += 1

    return paths


def generate_test_anchors() -> List[Anchor]:
    """Generate test anchor data with both forward and reverse strand hits."""
    anchors = []
    read_id = "test_read_001"

    # Forward strand hits
    # Path 0: Dense cluster in middle region
    for i in range(25):
        anchors.append(Anchor(
            read_id=read_id, path_id=0, path_name="HLA-DRB1*01:01",
            query_pos=500 + i * 30, ref_coord=7000 + i * 40, length=20
        ))

    # Path 1: Scattered anchors
    for i, ref in enumerate([1000, 2500, 4000, 5500, 7000]):
        anchors.append(Anchor(
            read_id=read_id, path_id=1, path_name="HLA-DRB1*03:01",
            query_pos=200 + i * 100, ref_coord=ref, length=20
        ))

    # Reverse strand hits
    # Path 5 (reverse of 01:01): Some hits
    for i in range(10):
        anchors.append(Anchor(
            read_id=read_id, path_id=5, path_name="HLA-DRB1*01:01_reverse",
            query_pos=100 + i * 50, ref_coord=3000 + i * 60, length=20
        ))

    # Path 6 (reverse of 03:01): Dense cluster
    for i in range(15):
        anchors.append(Anchor(
            read_id=read_id, path_id=6, path_name="HLA-DRB1*03:01_reverse",
            query_pos=300 + i * 25, ref_coord=10000 + i * 30, length=20
        ))

    return anchors


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_anchor_heatmap_grouped(anchors: List[Anchor], paths: List[Path], output_file: str,
                                 window_size: int = 200, stride: int = 50):
    """Generate heatmap with +/- sub-bars grouped by haplotype.

    Args:
        anchors: List of anchor mappings (single read)
        paths: List of haplotype paths (forward and reverse)
        output_file: Output file path for the plot
        window_size: Size of sliding window in bp (default: 200)
        stride: Stride/step size for sliding window in bp (default: 50)
    """
    # Group paths by haplotype
    groups = group_paths_by_haplotype(paths)

    if not groups:
        print("No paths to plot")
        return

    # Build path name to path mapping
    path_by_name = {p.name: p for p in paths}

    # Determine global coordinate range
    max_length = max(g.length for g in groups)
    num_windows = int(np.ceil((max_length - window_size) / stride)) + 1

    # Calculate figure height based on number of groups
    fig_height = max(6, len(groups) * 0.8 + 2)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    # Build density matrices for forward and reverse
    # We'll have 2 rows per haplotype group
    num_rows = len(groups) * 2
    density_matrix = np.ma.zeros((num_rows, num_windows))

    # Build mapping: row_idx -> (group, strand)
    row_info = []
    for group in groups:
        row_info.append((group, "+"))  # Forward
        row_info.append((group, "-"))  # Reverse

    # Populate density matrix
    for anchor in anchors:
        path = path_by_name.get(anchor.path_name)
        if not path:
            continue

        # Find row index for this path
        for row_idx, (group, strand) in enumerate(row_info):
            if path.strand == strand and path.base_name == group.base_name:
                # Anchor interval
                anchor_start = anchor.ref_coord
                anchor_end = anchor.ref_coord + anchor.length

                # Find overlapping windows
                first_window = max(0, int((anchor_start - window_size) / stride))
                last_window = min(num_windows - 1, int(anchor_end / stride))

                for win_idx in range(first_window, last_window + 1):
                    window_start = win_idx * stride
                    window_end = window_start + window_size
                    overlap_start = max(anchor_start, window_start)
                    overlap_end = min(anchor_end, window_end)
                    if overlap_start < overlap_end:
                        density_matrix[row_idx, win_idx] += 1
                break

    # Mask windows beyond each path's length
    for row_idx, (group, strand) in enumerate(row_info):
        path_length = group.length
        path_max_window = int(np.ceil((path_length - window_size) / stride))
        if path_max_window < num_windows - 1:
            density_matrix[row_idx, path_max_window + 1:] = np.ma.masked

    # Create colormap
    cmap = plt.cm.Blues
    cmap.set_bad(color='white')

    vmax = np.max(density_matrix) if np.max(density_matrix) > 0 else 1

    # Drawing parameters
    bar_height = 0.35  # Height of each +/- bar
    group_spacing = 0.3  # Space between haplotype groups

    y_positions = []  # Track y positions for labels
    y_pos = 0

    for group_idx, group in enumerate(groups):
        # Forward bar (+)
        forward_row = group_idx * 2
        y_forward = y_pos
        y_positions.append((y_forward, group, "+"))

        row_data = density_matrix[forward_row:forward_row + 1, :]
        im = ax.imshow(row_data,
                       aspect='auto',
                       cmap=cmap,
                       interpolation='nearest',
                       extent=[0, max_length, y_forward + bar_height, y_forward],
                       vmin=0, vmax=vmax, zorder=1)

        # Draw path length boundary
        ax.plot([group.length, group.length],
                [y_forward, y_forward + bar_height],
                color='dimgray', linewidth=1, linestyle='-', zorder=2)

        # Reverse bar (-)
        reverse_row = group_idx * 2 + 1
        y_reverse = y_forward + bar_height
        y_positions.append((y_reverse, group, "-"))

        row_data = density_matrix[reverse_row:reverse_row + 1, :]
        ax.imshow(row_data,
                  aspect='auto',
                  cmap=cmap,
                  interpolation='nearest',
                  extent=[0, max_length, y_reverse + bar_height, y_reverse],
                  vmin=0, vmax=vmax, zorder=1)

        ax.plot([group.length, group.length],
                [y_reverse, y_reverse + bar_height],
                color='dimgray', linewidth=1, linestyle='-', zorder=2)

        # Add haplotype name label (centered on group)
        # Split long names by ':' for multi-line display
        label_y = y_forward + bar_height  # Center between + and -
        if ':' in group.base_name:
            parts = group.base_name.split(':', 1)
            label_text = f"{parts[0]}\n{parts[1]}"
        else:
            label_text = group.base_name
        ax.text(-max_length * 0.01, label_y,
                label_text,
                va='center', ha='right', fontsize=9, fontweight='normal',
                linespacing=0.9)

        # Add +/- labels
        ax.text(-max_length * 0.005, y_forward + bar_height / 2,
                '+', va='center', ha='right', fontsize=8, color='green')
        ax.text(-max_length * 0.005, y_reverse + bar_height / 2,
                '-', va='center', ha='right', fontsize=8, color='red')

        # Add length annotation (once per group, on right side)
        ax.text(group.length + max_length * 0.01, label_y,
                f'{group.length:,} bp',
                va='center', ha='left', fontsize=9, color='dimgray',
                fontstyle='italic')

        # Move to next group with spacing
        y_pos = y_reverse + bar_height + group_spacing

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                        pad=0.08, fraction=0.04, aspect=30)
    cbar.set_label('Anchor Count per Bin', fontsize=10, fontweight='bold')

    # Formatting
    ax.set_xlabel("Reference Coordinate (bp)", fontsize=11)
    ax.set_ylabel("Haplotype Path", fontsize=11)

    # Add read ID to title
    read_id = anchors[0].read_id if anchors else "No anchors"
    ax.set_title(f"Anchor Density Heatmap: {read_id}\n(Window: {window_size} bp, Stride: {stride} bp)",
                 fontsize=12, pad=15)

    # Set axis limits
    ax.set_xlim(-max_length * 0.15, max_length * 1.12)
    ax.set_ylim(y_pos - group_spacing, -0.1)

    # Hide y-axis ticks (we use custom labels)
    ax.set_yticks([])

    # Clean background
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved heatmap to {output_file}")


def print_summary_statistics(anchors: List[Anchor], paths: List[Path]):
    """Print summary statistics about anchor distribution."""
    if not anchors:
        print("\nNo anchors to summarize.\n")
        return

    read_id = anchors[0].read_id
    total_anchors = len(anchors)

    print("\n" + "=" * 80)
    print(f"ANCHOR DISTRIBUTION SUMMARY: {read_id}")
    print("=" * 80)

    # Count anchors per path
    path_counts = defaultdict(int)
    for anchor in anchors:
        path_counts[anchor.path_name] += 1

    # Group by haplotype
    groups = group_paths_by_haplotype(paths)
    path_by_name = {p.name: p for p in paths}

    print(f"\nTotal anchors: {total_anchors}")
    print(f"Haplotypes: {len(groups)}")
    print(f"Paths hit: {len(path_counts)}/{len(paths)}\n")

    print("Per-Haplotype Distribution:")
    print("-" * 80)
    print(f"{'Haplotype':<40} {'Forward (+)':<15} {'Reverse (-)':<15} {'Total':<10}")
    print("-" * 80)

    for group in groups:
        fwd_count = 0
        rev_count = 0

        if group.forward_path:
            fwd_count = path_counts.get(group.forward_path.name, 0)
        if group.reverse_path:
            rev_count = path_counts.get(group.reverse_path.name, 0)

        total = fwd_count + rev_count
        pct = (total / total_anchors) * 100 if total_anchors > 0 else 0

        print(f"{group.base_name:<40} {fwd_count:<15} {rev_count:<15} {total:<10} ({pct:.1f}%)")

    print("-" * 80)
    print()


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Main entry point for the anchor heatmap tool."""
    parser = argparse.ArgumentParser(
        description="Visualize anchor density distribution across haplotype paths",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with built-in test data
  python plot_anchor_heatmap.py

  # Load from piru anchor dump
  python plot_anchor_heatmap.py --input anchors.tsv

  # Custom output and window settings
  python plot_anchor_heatmap.py --input anchors.tsv --output my_heatmap.png --window-size 300
"""
    )

    parser.add_argument(
        "--input", "-i",
        help="Input anchor dump file from piru (TSV format). If not provided, uses test data."
    )

    parser.add_argument(
        "--output", "-o",
        default="anchor_heatmap.png",
        help="Output file path for heatmap (default: anchor_heatmap.png)"
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=200,
        help="Sliding window size in bp (default: 200)"
    )

    parser.add_argument(
        "--stride",
        type=int,
        default=50,
        help="Sliding window stride in bp (default: 50)"
    )

    args = parser.parse_args()

    # Load data
    if args.input:
        print(f"Loading anchor dump from {args.input}...")
        paths, anchors = parse_anchor_dump(args.input)
    else:
        print("Generating test data...")
        paths = generate_test_paths()
        anchors = generate_test_anchors()

    read_id = anchors[0].read_id if anchors else "N/A"
    print(f"  Read ID: {read_id}")
    print(f"  Paths: {len(paths)}")
    print(f"  Anchors: {len(anchors)}")

    # Group and report
    groups = group_paths_by_haplotype(paths)
    print(f"  Haplotypes: {len(groups)}")

    # Plot heatmap
    print(f"\nGenerating grouped heatmap (window: {args.window_size} bp, stride: {args.stride} bp)...")
    plot_anchor_heatmap_grouped(anchors, paths, args.output,
                                 window_size=args.window_size, stride=args.stride)

    # Print summary statistics
    print_summary_statistics(anchors, paths)


if __name__ == "__main__":
    main()
