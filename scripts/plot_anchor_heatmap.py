#!/usr/bin/env python3
"""Plot heatmap of anchor positions across haplotype paths.

This tool visualizes where read anchors map across different haplotype paths,
helping diagnose whether anchors are properly distributed or concentrated on
a single path.

Supports grouped visualization of forward (+) and reverse (-) strands per haplotype.

Input can be either:
  - A single anchor dump file (legacy)
  - A directory containing per-read anchor files (<read_id>_anchors.tsv)
"""

import argparse
import os
import glob
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
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


@dataclass
class PafMapping:
    """Represents a single PAF alignment record."""
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str  # '+' or '-'
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    block_length: int
    mapq: int
    chain_score: Optional[int] = None  # From cs:i: or s1:i: tag


# =============================================================================
# Input Parsing
# =============================================================================

def parse_paf_file(filepath: str) -> Dict[str, List[PafMapping]]:
    """Parse a PAF file and return mappings grouped by query name.

    Args:
        filepath: Path to the PAF file

    Returns:
        Dictionary mapping query_name to list of PafMapping objects
    """
    mappings: Dict[str, List[PafMapping]] = defaultdict(list)

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 12:
                continue

            # Parse optional tags for chain score (cs:i:, s1:i:, or sc:i:)
            chain_score = None
            for tag in parts[12:]:
                if tag.startswith('cs:i:') or tag.startswith('s1:i:') or tag.startswith('sc:i:'):
                    try:
                        chain_score = int(tag.split(':')[2])
                        break
                    except (IndexError, ValueError):
                        pass

            mapping = PafMapping(
                query_name=parts[0],
                query_length=int(parts[1]),
                query_start=int(parts[2]),
                query_end=int(parts[3]),
                strand=parts[4],
                target_name=parts[5],
                target_length=int(parts[6]),
                target_start=int(parts[7]),
                target_end=int(parts[8]),
                num_matches=int(parts[9]),
                block_length=int(parts[10]),
                mapq=int(parts[11]),
                chain_score=chain_score
            )
            mappings[mapping.query_name].append(mapping)

    return mappings


def list_anchor_files(directory: str) -> List[Tuple[str, str]]:
    """List all anchor dump files in a directory.

    Returns:
        List of (read_num, filepath) tuples sorted numerically.
        Filename format: read_<N>_anchors.tsv
    """
    pattern = os.path.join(directory, "read_*_anchors.tsv")
    files = glob.glob(pattern)

    results = []
    for filepath in files:
        filename = os.path.basename(filepath)
        # Extract number from filename: read_<N>_anchors.tsv
        # Remove "read_" prefix and "_anchors.tsv" suffix
        num_str = filename[5:-12]  # len("read_") = 5, len("_anchors.tsv") = 12
        try:
            read_num = int(num_str)
            results.append((str(read_num), filepath))
        except ValueError:
            continue

    return sorted(results, key=lambda x: int(x[0]))


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
                                 window_size: int = 200, stride: int = 50, label: str = "",
                                 ax: Optional[plt.Axes] = None, fig: Optional[plt.Figure] = None,
                                 show_colorbar: bool = True, title_prefix: str = "",
                                 paf_mappings: Optional[List[PafMapping]] = None):
    """Generate heatmap with +/- sub-bars grouped by haplotype.

    Args:
        anchors: List of anchor mappings (single read)
        paths: List of haplotype paths (forward and reverse)
        output_file: Output file path for the plot (ignored if ax is provided)
        window_size: Size of sliding window in bp (default: 200)
        stride: Stride/step size for sliding window in bp (default: 50)
        label: Custom label text to display in top-left corner
        ax: Optional matplotlib Axes to plot on (for subplot mode)
        fig: Optional matplotlib Figure (for subplot mode)
        show_colorbar: Whether to show the colorbar (default: True)
        title_prefix: Optional prefix for the title (e.g., "A) " for subplots)
        paf_mappings: Optional list of PAF mappings for this read to draw as intervals
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

    # Determine if we're in standalone or subplot mode
    standalone_mode = (ax is None)

    if standalone_mode:
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
                # Anchor interval - convert coordinates for reverse strand
                if strand == "-":
                    # Flip coordinates for reverse complement: pos -> path_length - pos
                    path_length = group.length
                    anchor_start = path_length - (anchor.ref_coord + anchor.length)
                    anchor_end = path_length - anchor.ref_coord
                else:
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

    # Create custom blue colormap: light blue (0) → dark blue (30) → darker blue (max)
    from matplotlib.colors import LinearSegmentedColormap

    data_max = max(np.max(density_matrix), 40)  # Ensure at least 40

    color_stops = [
        (0 / data_max, '#f7fbff'),
        (7.5 / data_max, '#c6dbef'),
        (15 / data_max, '#6baed6'),
        (22.5 / data_max, '#2171b5'),
        (30 / data_max, '#08306b'),   # dark blue at 30
        (1.0, '#041529'),             # darker blue at max
    ]

    cmap = LinearSegmentedColormap.from_list(
        "anchor_blue",
        list(zip([c[0] for c in color_stops], [c[1] for c in color_stops]))
    )
    cmap.set_bad(color='white')

    vmax = data_max

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
                       vmin=0, vmax=vmax,
                       interpolation='nearest',
                       extent=[0, max_length, y_forward + bar_height, y_forward],
                       zorder=1)

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
                  vmin=0, vmax=vmax,
                  interpolation='nearest',
                  extent=[0, max_length, y_reverse + bar_height, y_reverse],
                  zorder=1)

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

    # Draw PAF mapping intervals if provided
    if paf_mappings:
        # Build lookup from (base_name, strand) -> y_position
        path_to_y: Dict[Tuple[str, str], float] = {}
        for y, group, strand in y_positions:
            path_to_y[(group.base_name, strand)] = y

        # Use PAF order (already sorted by chain score rank)
        # Draw each mapping as an interval line
        for idx, mapping in enumerate(paf_mappings, start=1):
            # Use target_name as base_name and strand from PAF field
            base_name = mapping.target_name
            strand = mapping.strand  # '+' or '-' from PAF column 5

            # Get the y position and group for this path
            y = path_to_y.get((base_name, strand))
            if y is None:
                # Try matching without strand (in case PAF strand differs)
                continue

            # PAF coordinates are always in forward reference space, no flipping needed
            display_start = mapping.target_start
            display_end = mapping.target_end

            # Draw interval line - centered in the middle of the bar
            interval_y = y + bar_height / 2

            line_color = '#d97706'
            line_width = 1.5
            line_halo = [pe.withStroke(linewidth=3.5, foreground='white', alpha=0.85)]

            # Draw the interval: |------|
            ax.plot([display_start, display_end], [interval_y, interval_y],
                    color=line_color, linewidth=line_width, solid_capstyle='butt', zorder=10,
                    path_effects=line_halo)
            # Left cap |
            ax.plot([display_start, display_start],
                    [interval_y - 0.03, interval_y + 0.03],
                    color=line_color, linewidth=line_width, zorder=10,
                    path_effects=line_halo)
            # Right cap |
            ax.plot([display_end, display_end],
                    [interval_y - 0.03, interval_y + 0.03],
                    color=line_color, linewidth=line_width, zorder=10,
                    path_effects=line_halo)

            # Add label inside the chain interval, centered with rounded box
            label_x = (display_start + display_end) / 2
            if mapping.chain_score is not None:
                label_text = f"C{idx} {mapping.chain_score}"
            else:
                label_text = f"C{idx}"
            ax.text(label_x, interval_y, label_text,
                    ha='center', va='center', fontsize=7, fontweight='medium',
                    color='#92400e', zorder=12,
                    bbox=dict(boxstyle='round,pad=0.2,rounding_size=0.3',
                              facecolor='white', edgecolor='none'))

    # Add colorbar (only in standalone mode or if explicitly requested)
    if show_colorbar and standalone_mode:
        cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                            pad=0.08, fraction=0.04, aspect=30)
        cbar.set_label('Anchor Count per Bin', fontsize=10, fontweight='bold')

    # Formatting
    ax.set_xlabel("Reference Coordinate (bp)", fontsize=11 if standalone_mode else 9)
    ax.set_ylabel("Haplotype Path", fontsize=11 if standalone_mode else 9)

    # Add read ID to title
    read_id = anchors[0].read_id if anchors else "No anchors"
    if standalone_mode:
        ax.set_title(f"Anchor Density Heatmap: {read_id}\n(Window: {window_size} bp, Stride: {stride} bp)",
                     fontsize=12, pad=15)
    else:
        # Shorter title for subplot mode
        ax.set_title(f"{title_prefix}{read_id}", fontsize=10, pad=8)

    # Set axis limits with padding for interval labels
    ax.set_xlim(-max_length * 0.15, max_length * 1.12)
    ax.set_ylim(y_pos - group_spacing + 0.5, -0.6)

    # Hide y-axis ticks (we use custom labels)
    ax.set_yticks([])

    # Clean background
    ax.set_facecolor('white')
    if fig is not None:
        fig.patch.set_facecolor('white')

    # Add custom label if provided (top-left of figure, outside plot)
    if label and standalone_mode and fig is not None:
        fig.text(0.01, 0.99, label, fontsize=11, fontstyle='italic',
                 color='goldenrod', va='top', ha='left')

    # Only save if in standalone mode
    if standalone_mode:
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
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


def plot_combined_heatmaps(anchor_files: List[Tuple[str, str]], output_file: str,
                           window_size: int, stride: int, label: str = "",
                           paf_mappings_by_read: Optional[Dict[str, List[PafMapping]]] = None):
    """Plot all reads as subplots in a single figure, sorted by read name.

    Args:
        anchor_files: List of (read_num, filepath) tuples
        output_file: Output file path for the combined plot
        window_size: Sliding window size in bp
        stride: Sliding window stride in bp
        label: Custom label text to display
        paf_mappings_by_read: Optional dict mapping read_id to list of PAF mappings
    """
    # First pass: parse all files to get read_id for sorting
    read_data = []  # List of (read_id, paths, anchors, filepath)

    print("Parsing all anchor files...")
    for read_num, filepath in anchor_files:
        paths, anchors = parse_anchor_dump(filepath)
        if not anchors:
            print(f"  Warning: No anchors in {filepath}, skipping")
            continue
        read_id = anchors[0].read_id
        read_data.append((read_id, paths, anchors, filepath))
        print(f"  {filepath}: {read_id} ({len(anchors)} anchors)")

    if not read_data:
        print("Error: No valid anchor files found")
        return False

    # Sort by read_id (alphabetically)
    read_data.sort(key=lambda x: x[0])

    print(f"\nSorted {len(read_data)} reads by read name")

    # Determine subplot layout: one read per row for easy scrolling
    n_reads = len(read_data)
    n_cols = 1
    n_rows = n_reads

    # Calculate figure size
    fig_width = 16
    # Get number of haplotypes from first read to estimate height
    groups = group_paths_by_haplotype(read_data[0][1])
    subplot_height = max(4, len(groups) * 0.6 + 1.5)
    fig_height = subplot_height * n_rows + 1  # Extra space for suptitle

    print(f"Creating figure with {n_rows}x{n_cols} subplots...")
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))

    # Flatten axes array for easier iteration
    if n_reads == 1:
        axes = [axes]
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    else:
        axes = axes.flatten()

    # Plot each read
    for idx, (read_id, paths, anchors, filepath) in enumerate(read_data):
        ax = axes[idx]
        # Get PAF mappings for this read if available
        read_paf_mappings = None
        if paf_mappings_by_read:
            read_paf_mappings = paf_mappings_by_read.get(read_id)
        print(f"  Plotting [{idx+1}/{n_reads}] {read_id}..." +
              (f" ({len(read_paf_mappings)} PAF mappings)" if read_paf_mappings else ""))
        plot_anchor_heatmap_grouped(
            anchors, paths, output_file,
            window_size=window_size, stride=stride,
            ax=ax, fig=fig, show_colorbar=False,
            title_prefix="",
            paf_mappings=read_paf_mappings
        )

    # Hide empty subplots
    for idx in range(n_reads, len(axes)):
        axes[idx].set_visible(False)

    # Add overall title
    fig.suptitle(f"Anchor Density Heatmaps (Window: {window_size} bp, Stride: {stride} bp)",
                 fontsize=14, fontweight='bold', y=0.99)

    # Add custom label if provided
    if label:
        fig.text(0.01, 0.99, label, fontsize=11, fontstyle='italic',
                 color='goldenrod', va='top', ha='left')

    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Leave room for suptitle
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"\nSaved combined heatmap to {output_file}")
    return True


# =============================================================================
# Main Entry Point
# =============================================================================

def process_single_file(filepath: str, output_file: str, window_size: int, stride: int,
                         label: str = "",
                         paf_mappings_by_read: Optional[Dict[str, List[PafMapping]]] = None):
    """Process a single anchor dump file and generate heatmap."""
    print(f"Loading anchor dump from {filepath}...")
    paths, anchors = parse_anchor_dump(filepath)

    if not anchors:
        print(f"  Warning: No anchors found, skipping")
        return False

    read_id = anchors[0].read_id
    print(f"  Read ID: {read_id}")
    print(f"  Paths: {len(paths)}")
    print(f"  Anchors: {len(anchors)}")

    # Look up PAF mappings for this read
    paf_mappings = None
    if paf_mappings_by_read:
        paf_mappings = paf_mappings_by_read.get(read_id)
        if paf_mappings:
            print(f"  PAF mappings: {len(paf_mappings)}")

    # Group and report
    groups = group_paths_by_haplotype(paths)
    print(f"  Haplotypes: {len(groups)}")

    # Plot heatmap
    print(f"  Generating heatmap (window: {window_size} bp, stride: {stride} bp)...")
    plot_anchor_heatmap_grouped(anchors, paths, output_file,
                                 window_size=window_size, stride=stride, label=label,
                                 paf_mappings=paf_mappings)

    # Print summary statistics
    print_summary_statistics(anchors, paths)
    return True


def main():
    """Main entry point for the anchor heatmap tool."""
    parser = argparse.ArgumentParser(
        description="Visualize anchor density distribution across haplotype paths",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with built-in test data
  python plot_anchor_heatmap.py

  # Load from single anchor dump file
  python plot_anchor_heatmap.py --input read_0_anchors.tsv

  # Process all files in a directory (generates one heatmap per read)
  python plot_anchor_heatmap.py --input /tmp/anchors

  # Combine all reads into a single image (sorted by read name)
  python plot_anchor_heatmap.py --input /tmp/anchors --combined -o combined.png

  # Custom window settings
  python plot_anchor_heatmap.py --input /tmp/anchors --window-size 300

  # Add PAF mapping intervals to heatmap
  python plot_anchor_heatmap.py --input read_0_anchors.tsv -p mappings.paf

  # Combined view with PAF intervals
  python plot_anchor_heatmap.py --input /tmp/anchors --combined -p mappings.paf
"""
    )

    parser.add_argument(
        "--input", "-i",
        help="Input anchor dump file or directory from piru. If not provided, uses test data."
    )

    parser.add_argument(
        "--output", "-o",
        default="anchor_heatmap.png",
        help="Output file path for heatmap (default: anchor_heatmap.png). "
             "For directory input, this is ignored and files are named read_<N>_heatmap.png"
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

    parser.add_argument(
        "--label", "-l",
        default="",
        help="Custom label text to display in top-left corner of the plot"
    )

    parser.add_argument(
        "--combined", "-c",
        action="store_true",
        help="Combine all reads into a single image as subplots (sorted by read name). "
             "Only applies when --input is a directory."
    )

    parser.add_argument(
        "--paf", "-p",
        help="Optional PAF file with mapping results. When provided, draws interval markers "
             "on the heatmap showing mapping locations. Multiple mappings are numbered (1, 2, 3...)."
    )

    args = parser.parse_args()

    # Parse PAF file if provided
    paf_mappings_by_read: Dict[str, List[PafMapping]] = {}
    if args.paf:
        print(f"Loading PAF mappings from {args.paf}...")
        paf_mappings_by_read = parse_paf_file(args.paf)
        total_mappings = sum(len(m) for m in paf_mappings_by_read.values())
        print(f"  Loaded {total_mappings} mappings for {len(paf_mappings_by_read)} reads")

    # Load data
    if args.input:
        if os.path.isdir(args.input):
            # Directory mode - process all files
            anchor_files = list_anchor_files(args.input)

            if not anchor_files:
                print(f"Error: No anchor files found in {args.input}")
                print("Expected files matching pattern: read_*_anchors.tsv")
                return 1

            print(f"Found {len(anchor_files)} anchor files in {args.input}")
            print("=" * 60)

            if args.combined:
                # Combined mode: all reads in one image, sorted by read name
                # Output in input directory, named after directory unless explicitly specified
                if args.output != "anchor_heatmap.png":
                    output_file = args.output
                else:
                    dir_name = os.path.basename(os.path.normpath(args.input))
                    output_file = os.path.join(args.input, f"{dir_name}.png")
                success = plot_combined_heatmaps(
                    anchor_files, output_file,
                    args.window_size, args.stride, args.label,
                    paf_mappings_by_read=paf_mappings_by_read if paf_mappings_by_read else None
                )
                return 0 if success else 1
            else:
                # Individual mode: one heatmap per read
                success_count = 0
                for read_num, filepath in anchor_files:
                    output_file = os.path.join(args.input, f"read_{read_num}_heatmap.png")
                    print(f"\n[{int(read_num)+1}/{len(anchor_files)}] Processing read_{read_num}...")
                    if process_single_file(filepath, output_file, args.window_size, args.stride, args.label,
                                           paf_mappings_by_read=paf_mappings_by_read if paf_mappings_by_read else None):
                        success_count += 1

                print("\n" + "=" * 60)
                print(f"Done. Generated {success_count}/{len(anchor_files)} heatmaps.")
                return 0

        else:
            # Single file mode - name output after input file unless explicitly specified
            if args.output != "anchor_heatmap.png":
                output_file = args.output
            else:
                output_file = f"{args.input}.png"
            return 0 if process_single_file(args.input, output_file, args.window_size, args.stride, args.label,
                                            paf_mappings_by_read=paf_mappings_by_read if paf_mappings_by_read else None) else 1

    else:
        # Test data mode
        print("Generating test data...")
        paths = generate_test_paths()
        anchors = generate_test_anchors()

        read_id = anchors[0].read_id if anchors else "N/A"
        print(f"  Read ID: {read_id}")
        print(f"  Paths: {len(paths)}")
        print(f"  Anchors: {len(anchors)}")

        groups = group_paths_by_haplotype(paths)
        print(f"  Haplotypes: {len(groups)}")

        print(f"\nGenerating grouped heatmap (window: {args.window_size} bp, stride: {args.stride} bp)...")
        plot_anchor_heatmap_grouped(anchors, paths, args.output,
                                     window_size=args.window_size, stride=args.stride, label=args.label)
        print_summary_statistics(anchors, paths)
        return 0


if __name__ == "__main__":
    import sys
    sys.exit(main() or 0)
