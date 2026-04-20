#!/usr/bin/env python3
"""
Generate annotated tree figure for metagenome bins.
Automatically selects circular layout for >50 bins.
Designed to be called from metagenome_context.py

python scripts/metagenome_tree_figure.py \
    --csv '/Users/kh18/server_data/metagenomes/kaTriNubi1/bin_data.csv' \
    --outdir /Users/kh18/Documents/md_ASG/Trididemnum_nubilum \
    --filename metagenome_tree.png 

    
python scripts/metagenome_taxonomy_tree.py \
    --csv  /Users/kh18/server_data/metagenomes/ohBolCyan1/bin_data.csv \
    --outdir Bolosoma_cyanae/tree_assets \
    --use-taxid --ncbi-db ~/ncbi_taxadb.sqlite


python scripts/metagenome_tree_figure.py \
    --csv /Users/kh18/server_data/metagenomes/ohBolCyan1/bin_data.csv \
    --outdir Bolosoma_cyanae/ \
    --filename metagenome_tree.png \
    --order-json Bolosoma_cyanae/tree_assets/taxonomy_tree.json

"""

import argparse
import json
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, Rectangle, Wedge
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager
import warnings
warnings.filterwarnings('ignore')

try:
    from ete3 import Tree as EteTree
except Exception:
    EteTree = None


CIRCULAR_SPAN_FRACTION = 350.0 / 360.0


def _spread_angles(n_points, span_fraction=CIRCULAR_SPAN_FRACTION):
    """Return evenly spaced angles that leave a gap in the circle."""
    if n_points <= 0:
        return np.array([]), 0.0, 0.0, 0.0
    span_fraction = min(max(span_fraction, 0.1), 1.0)
    angle_span = 2 * np.pi * span_fraction
    gap_angle = 2 * np.pi - angle_span
    offset = gap_angle / 2
    angles = offset + np.linspace(0, angle_span, n_points, endpoint=False)
    angle_step = angle_span / n_points
    return angles, angle_step, offset, angle_span


def _interpolate_angles(start, end, num=30):
    """Return angle samples along the shortest arc between start and end."""
    delta = end - start
    if delta > np.pi:
        delta -= 2 * np.pi
    elif delta < -np.pi:
        delta += 2 * np.pi
    return start + np.linspace(0, delta, num)


def resolve_open_sans_font(env_var="GENOMENOTES_FONT"):
    """Locate a regular Open Sans font file."""
    import os
    
    # 1. explicit override
    p = os.environ.get(env_var)
    if p and Path(p).is_file():
        return p
    
    # helper to pick first non-italic from list of paths
    def pick_upright(paths):
        regular = [x for x in paths if "Regular" in str(x)]
        if regular:
            return str(sorted(regular)[0])
        upright = [x for x in paths if "italic" not in x.name.lower()]
        if upright:
            return str(sorted(upright)[0])
        return str(sorted(paths)[0]) if paths else None
    
    # Check common font directories
    font_dirs = [
        Path.home() / "Library" / "Fonts",  # macOS
        Path("/usr/share/fonts/truetype/open-sans"),
        Path("/usr/share/fonts/truetype"),
        Path("/usr/local/share/fonts"),
        Path.home() / ".local/share/fonts",
        Path.home() / ".fonts"
    ]
    
    for font_dir in font_dirs:
        if font_dir.exists():
            hits = []
            for pattern in ["OpenSans-Regular.ttf", "OpenSans*.ttf", "open-sans*.ttf"]:
                hits.extend(font_dir.glob(pattern))
            chosen = pick_upright(hits)
            if chosen:
                return chosen
    
    return None


def setup_font():
    """Configure OpenSans font if available."""
    try:
        font_path = resolve_open_sans_font()
        if font_path:
            font_manager.fontManager.addfont(font_path)
            plt.rcParams['font.family'] = 'Open Sans'
            plt.rcParams['font.sans-serif'] = ['Open Sans', 'DejaVu Sans', 'sans-serif']
            return True
    except Exception:
        pass
    return False


def format_species_name(name):
    r"""Format taxonomic names in italics.
    
    Uses matplotlib's text styling to italicize without collapsing spaces.
    """
    if not isinstance(name, str):
        return name
    
    # Use style parameter instead - simpler and preserves spaces
    # But we need to return the name and apply style='italic' in the text() call
    # OR we can use the fontproperties approach
    # OR split and wrap each word separately
    
    # Split by spaces and wrap each word, preserving spaces
    words = name.split()
    formatted_words = [rf"$\mathit{{{word}}}$" for word in words]
    return ' '.join(formatted_words)


def load_order_metadata(order_json):
    """Read bin ordering and optional tree path from taxonomy_tree.json."""
    path = Path(order_json)
    if not path.is_file():
        raise FileNotFoundError(f"Order file not found: {order_json}")
    data = json.loads(path.read_text())
    order = data.get('bin_order')
    if not isinstance(order, list) or not order:
        raise ValueError(f"No 'bin_order' list found in {order_json}")
    tree_path = data.get('tree_newick')
    if tree_path:
        candidate = Path(tree_path)
        if candidate.is_absolute():
            tree_path = candidate
        elif candidate.exists():
            tree_path = candidate.resolve()
        else:
            tree_path = (path.parent / candidate).resolve()
    return order, tree_path


class MetagenomeTreeFigure:
    """Generate phylogenetic tree figures for metagenome bins."""

    
    def __init__(self, csv_path, bin_order=None, tree_path=None):
        self.csv_path = csv_path
        self.df = None
        self.colors = self._setup_colors()
        self.bin_order = bin_order or []
        self.tree_path = tree_path
        self.tree = None
        self.tree_layout = None
        setup_font()
        
    
    def _setup_colors(self):
        """Setup color schemes."""
        merianbow4 = [
            "#710093", "#0080ff","#009286","#A676FF","#EE6A15", "#BC007B",
            "#ABA9E5", "#EA005B", "#E9CB19", "#fc1d1d", "#00C8FF",
            "#FD9514", "#87BF13", "#00AB3E" ,"#A1D0FF", "#006DDB", 
            "#ff1aef", "#FF6B70", "#A16C00", "#4E6400",  "#00e251", 
            "#00B0E0", "#AE83B6", "#005200", "#666666",
            "#A56183", "#8E454F", "#6B3122"
        ]
        
        # Red gradient for coverage; blue for completeness
        coverage_cmap = LinearSegmentedColormap.from_list(
            'cov', ['#fff5f5', '#fca5a5', '#b91c1c'], N=100
        )
        completeness_cmap = LinearSegmentedColormap.from_list(
            'comp', ['#f7fbff', '#6baed6', '#08519c'], N=100
        )
        
        return {
            'phylum_palette': merianbow4,
            'mag_color': '#c7c7c7',           # light grey for MAGs
            'circular_mag_color': '#000000',  # black for fully circular MAGs
            'coverage_cmap': coverage_cmap,
            'completeness_cmap': completeness_cmap,
            'size_color': '#888888'
        }

    def load_and_process(self):
        """Load and process the bin data."""
        self.df = pd.read_csv(self.csv_path)
        
        # Clean bin_id
        if 'bin_id' in self.df.columns:
            self.df['bin_id'] = self.df['bin_id'].astype(str).str.strip()
            self.df = self.df.drop_duplicates(subset='bin_id', keep='first')
        else:
            # still drop exact duplicate rows if no explicit bin identifier
            self.df = self.df.drop_duplicates()
        
        # Convert numeric columns
        numeric_cols = ['size', 'contigs', 'circular', 'mean_coverage', 
                       'unique_trnas', 'Contamination', 'Completeness']
        for col in numeric_cols:
            if col in self.df.columns:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
        
        # Clean rRNA columns
        for col in ['rrna_5s', 'rrna_16s', 'rrna_23s']:
            if col in self.df.columns:
                self.df[col] = self.df[col].astype(str).str.upper().str.strip()
        
        # Parse taxonomy
        self._parse_taxonomy()
        
        # Identify MAGs
        self._identify_mags()
        
        # Create labels
        self._create_labels()
        self._maybe_enforce_external_order()
        self._load_tree()

        return len(self.df)

    def _maybe_enforce_external_order(self):
        """Apply an externally supplied bin order if provided."""
        if not self.bin_order or 'bin_id' not in self.df.columns:
            return
        present = set(self.df['bin_id'])
        filtered = [bid for bid in self.bin_order if bid in present]
        missing = [bid for bid in self.df['bin_id'] if bid not in filtered]
        applied_order = filtered + missing
        order_lookup = {bid: idx for idx, bid in enumerate(applied_order)}
        fallback = len(order_lookup)
        self.df['__external_order'] = self.df['bin_id'].map(lambda x: order_lookup.get(x, fallback))
        self.df.sort_values('__external_order', kind='mergesort', inplace=True)
        self.df.drop(columns='__external_order', inplace=True)
        self.bin_order = applied_order

    def _build_rectangular_tree_layout(self, branch_end, y_spacing=1.0):
        """Return Cartesian coordinates for the ete3 tree to drive the rectangular layout."""
        if not self.tree or 'bin_id' not in self.df.columns:
            return None

        tree = self.tree.copy()
        leaves = list(tree.iter_leaves())
        if not leaves:
            return None

        df_order = list(self.df['bin_id'])
        if df_order:
            order_lookup = {bid: idx for idx, bid in enumerate(df_order)}
            leaves.sort(key=lambda leaf: order_lookup.get(leaf.name, len(df_order)))
        else:
            tree.ladderize()
            leaves = list(tree.iter_leaves())

        leaf_order = [leaf.name for leaf in leaves]
        leaf_y = {name: idx * y_spacing for idx, name in enumerate(leaf_order)}

        max_depth = max(tree.get_distance(leaf) for leaf in leaves)
        if max_depth <= 0:
            max_depth = 1.0
        depth_scale = branch_end / max_depth

        coords = {}

        def assign(node):
            """Recursively assign coordinates to internal and leaf nodes."""
            if node.is_leaf():
                y_val = leaf_y.get(node.name, 0.0)
            else:
                child_coords = [assign(child) for child in node.children if child]
                y_val = float(np.mean([pt[1] for pt in child_coords])) if child_coords else 0.0
            depth = node.get_distance(tree)
            x_val = branch_end if node.is_leaf() else depth * depth_scale
            coords[node] = (x_val, y_val)
            return coords[node]

        assign(tree)

        edges = []
        for node in tree.traverse():
            for child in node.children:
                px, py = coords[node]
                cx, cy = coords[child]
                edges.append({'px': px, 'py': py, 'cx': cx, 'cy': cy})

        return {
            'leaf_order': leaf_order,
            'leaf_y': leaf_y,
            'edges': edges,
            'y_min': min(leaf_y.values()),
            'y_max': max(leaf_y.values()),
            'branch_end': branch_end,
        }

    def _load_tree(self):
        """Load the ete3 tree if a path was provided."""
        if not self.tree_path or EteTree is None:
            return
        try:
            tree = EteTree(str(self.tree_path), format=1)
        except Exception as exc:
            print(f"[metagenome_tree_figure] Failed to load tree {self.tree_path}: {exc}")
            return
        present = set(self.df['bin_id'])
        for leaf in list(tree.iter_leaves()):
            if leaf.name not in present:
                leaf.detach()
        if not list(tree.iter_leaves()):
            print(f"[metagenome_tree_figure] Tree {self.tree_path} has no overlapping leaves")
            return
        tree.ladderize()
        self.tree = tree

    def _build_tree_layout(self, inner_radius, tree_outer_radius):
        """Compute polar coordinates for the ete3 tree if available."""
        if not self.tree:
            return None
        tree = self.tree.copy()
        leaves = list(tree.iter_leaves())
        if not leaves:
            return None

        angles, _, _, _ = _spread_angles(len(leaves))
        leaf_angles = {}
        leaf_order = []
        for idx, leaf in enumerate(leaves):
            angle = angles[idx]
            leaf_angles[leaf.name] = angle
            leaf_order.append(leaf.name)

        max_depth = max(leaf.get_distance(tree) for leaf in leaves)
        if max_depth == 0:
            max_depth = 1.0
        radial_scale = (tree_outer_radius - inner_radius) / max_depth

        node_angles = {}
        node_radii = {}

        def assign(node):
            if node.is_leaf():
                angle = leaf_angles.get(node.name, 0.0)
            else:
                child_angles = [assign(child) for child in node.children if child]
                angle = float(np.mean(child_angles)) if child_angles else 0.0
            node_angles[node] = angle
            node_radii[node] = inner_radius + node.get_distance(tree) * radial_scale
            return angle

        assign(tree)

        edges = []
        for node in tree.traverse():
            if node.up is None:
                continue
            parent = node.up
            edges.append({
                'parent_angle': node_angles[parent],
                'parent_radius': node_radii[parent],
                'child_angle': node_angles[node],
                'child_radius': node_radii[node],
            })

        leaf_radii = {leaf.name: node_radii[leaf] for leaf in leaves}

        return {
            'leaf_angles': leaf_angles,
            'leaf_radii': leaf_radii,
            'edges': edges,
            'leaf_order': leaf_order,
            'tree_outer_radius': tree_outer_radius,
        }
    
    def _parse_taxonomy(self):
        """Parse taxonomy strings."""
        tax_col = None
        for col in ['ncbi_classification', 'classification']:
            if col in self.df.columns:
                tax_col = col
                break
        
        if tax_col:
            self.df['phylum'] = self.df[tax_col].str.extract(r'p__([^;]+)')
            self.df['phylum'] = self.df['phylum'].fillna('Unclassified')
            
            # Extract other levels
            levels = {
                'd__': 'domain',
                'c__': 'class',
                'o__': 'order',
                'f__': 'family',
                'g__': 'genus',
                's__': 'species'
            }
            
            for prefix, level in levels.items():
                self.df[level] = self.df[tax_col].apply(
                    lambda x: self._extract_tax_level(x, prefix) if pd.notna(x) else 'Unknown'
                )
    
    def _extract_tax_level(self, tax_string, prefix):
        """Extract specific taxonomic level from string."""
        if pd.isna(tax_string):
            return 'Unknown'
        parts = tax_string.split(';')
        for part in parts:
            if part.strip().startswith(prefix):
                value = part.strip()[3:]
                return value if value else 'Unknown'
        return 'Unknown'
    
    def _identify_mags(self):
        """Identify MAGs based on standard criteria."""
        def is_mag(row):
            bt = str(row.get('bin_type', '')).strip().lower()
            if bt:
                return 'mag' in bt

            has_rrnas = (
                row.get('rrna_5s', 'N') == 'Y' and
                row.get('rrna_16s', 'N') == 'Y' and
                row.get('rrna_23s', 'N') == 'Y'
            )
            try:
                trna_count = float(row.get('unique_trnas', 0))
            except Exception:
                trna_count = 0
            has_trnas = trna_count >= 18

            try:
                contamination = float(row.get('Contamination', 100))
            except Exception:
                contamination = 100
            contamination_ok = contamination <= 5

            try:
                completeness = float(row.get('Completeness', 0))
            except Exception:
                completeness = 0
            completeness_high = completeness >= 90

            # Require fully circular chromosomes for mid-completeness MAGs.
            try:
                contigs = float(row.get('contigs', 0))
                circular = float(row.get('circular', 0))
                fully_circular = (
                    (contigs == 1 and circular == 1) or
                    (contigs > 1 and contigs == circular)
                )
            except Exception:
                fully_circular = False

            completeness_mid_and_circular = completeness >= 50 and fully_circular
            
            return (
                contamination_ok and has_rrnas and has_trnas and
                (completeness_high or completeness_mid_and_circular)
            )
        
        self.df['is_MAG'] = self.df.apply(is_mag, axis=1)
        
        def is_fully_circular(row):
            if not row.get('is_MAG', False):
                return False
            contigs = row.get('contigs', 0)
            circular = row.get('circular', 0)
            return ((contigs == 1 and circular == 1) or 
                   (contigs > 1 and contigs == circular))
        
        self.df['is_fully_circular'] = self.df.apply(is_fully_circular, axis=1)
    
    def _create_labels(self):
        """Create informative labels from taxonomy."""
        if 'ncbi_taxon' in self.df.columns:
            self.df['label'] = self.df['ncbi_taxon']
        else:
            def create_label(row):
                for level in ['species', 'genus', 'family', 'order', 'class', 'phylum']:
                    if level in row and row[level] not in ['Unknown', '', None]:
                        label = row[level]
                        if level != 'species':
                            if 'domain' in row and row['domain'] == 'Archaea':
                                label += ' archaeon'
                            else:
                                label += ' bacterium'
                        return label
                return row.get('bin_id', 'Unknown')
            
            self.df['label'] = self.df.apply(create_label, axis=1)
    
    def create_figure(self, output_path, dpi=300):
        """Create the appropriate figure based on bin count."""
        n_bins = len(self.df)
        # Here is the threshold for circular vs rectangular
        if n_bins > 50:
            return self._create_circular_tree(output_path, dpi)
        else:
            print(f"[MetagenomeTreeFigure] rectangular layout selected for {n_bins} bins")
            return self._create_rectangular_tree(output_path, dpi)
    
    def _create_circular_tree(self, output_path, dpi=300):
        """Create circular tree for >50 bins."""
        fig, ax = plt.subplots(figsize=(14, 14), subplot_kw=dict(projection='polar'))

        # Sort by taxonomy only when no explicit order was provided
        if not self.bin_order:
            for col in ['phylum', 'class', 'order', 'family']:
                if col not in self.df.columns:
                    self.df[col] = 'Unclassified' if col == 'phylum' else 'Unknown'
            self.df = self.df.sort_values(['phylum', 'class', 'order', 'family'], na_position='last')

        # Setup colors per phylum
        phylum_colors = {}
        phylum_list = self.df['phylum'].unique()
        for i, phylum in enumerate(phylum_list):
            phylum_colors[phylum] = self.colors['phylum_palette'][i % len(self.colors['phylum_palette'])]

        # Calculate angles and radial layout constants
        n_bins = len(self.df)
        if n_bins == 0:
            raise ValueError("No bins available to plot")
        default_angles, angle_step, angle_offset, angle_span = _spread_angles(n_bins)

        inner_radius = 0.3
        branch_end_radius = 0.72
        tree_outer_radius = branch_end_radius - 0.12
        if tree_outer_radius <= inner_radius:
            tree_outer_radius = inner_radius + (branch_end_radius - inner_radius) * 0.6
        mag_radius = branch_end_radius
        comp_radius = branch_end_radius + 0.05
        cov_radius = branch_end_radius + 0.10
        size_base_radius = branch_end_radius + 0.16
        size_track_height = 0.18
        phylum_ring_inner = size_base_radius + size_track_height + 0.04
        phylum_ring_width = 0.05
        phylum_label_radius = phylum_ring_inner + phylum_ring_width + 0.03
        label_radius = size_base_radius + size_track_height + 0.02
        plot_max_radius = phylum_label_radius + 0.05
        show_tip_labels = n_bins <= 30

        tree_layout = self._build_tree_layout(inner_radius, tree_outer_radius)
        self.tree_layout = tree_layout
        if tree_layout:
            angles = []
            angle_lookup = tree_layout['leaf_angles']
            for idx, bin_id in enumerate(self.df['bin_id']):
                angle = angle_lookup.get(bin_id)
                if angle is None:
                    angle = default_angles[idx]
                angles.append(angle)
            angles = np.array(angles)
        else:
            angles = default_angles

        # Genome size scaling (Mbp)
        if 'size' in self.df.columns:
            max_size_val = pd.to_numeric(self.df['size'], errors='coerce').max()
            max_size_mbp = (max_size_val / 1e6) if pd.notna(max_size_val) else 0
        else:
            max_size_mbp = 0
        size_track_enabled = max_size_mbp > 0
        size_bar_width = angle_step * 0.6

        # Track phylum segments for the outer color ring
        phylum_segments = []
        current_phylum = None
        phylum_start_idx = 0

        if tree_layout:
            for edge in tree_layout['edges']:
                pa = edge['parent_angle']
                pr = edge['parent_radius']
                ca = edge['child_angle']
                cr = edge['child_radius']
                ax.plot([pa, pa], [pr, cr], color='#b5b5b5', linewidth=0.8, alpha=0.9)
                angle_vals = _interpolate_angles(pa, ca, num=32)
                ax.plot(angle_vals, np.full_like(angle_vals, cr), color='#b5b5b5', linewidth=0.8, alpha=0.9)

        for idx, (_, row) in enumerate(self.df.iterrows()):
            angle = angles[idx]

            if current_phylum is None:
                current_phylum = row['phylum']
                phylum_start_idx = idx
            elif row['phylum'] != current_phylum:
                phylum_segments.append((current_phylum, phylum_start_idx, idx))
                current_phylum = row['phylum']
                phylum_start_idx = idx

            # Draw branch
            if tree_layout:
                leaf_radius = tree_layout['leaf_radii'].get(row['bin_id'], tree_outer_radius)
                ax.plot([angle, angle], [leaf_radius, branch_end_radius],
                        color='#b0b0b0', linewidth=0.5, alpha=0.6)
            else:
                ax.plot([angle, angle], [inner_radius, branch_end_radius],
                        'k-', linewidth=0.5, alpha=0.3)

            # MAG indicator
            if row.get('is_MAG', False):
                if row.get('is_fully_circular', False):
                    ax.scatter(angle, mag_radius, s=80, c=[self.colors['circular_mag_color']],
                               marker='o', edgecolors='none', zorder=5)
                else:
                    ax.scatter(angle, mag_radius, s=80, c=[self.colors['mag_color']],
                               marker='o', edgecolors='none', zorder=5)

            # Quality indicators
            completeness = row.get('Completeness', 0)
            comp_color = self.colors['completeness_cmap'](max(min(completeness, 100), 0) / 100)
            ax.scatter(angle, comp_radius, s=60, c=[comp_color],
                       marker='s', edgecolors='black', linewidth=0.5)

            coverage = row.get('mean_coverage', 1)
            if coverage > 0:
                log_cov = np.log10(coverage)
                cov_norm = min(max(log_cov / 3, 0), 1)
            else:
                cov_norm = 0
            cov_color = self.colors['coverage_cmap'](cov_norm)
            ax.scatter(angle, cov_radius, s=60, c=[cov_color],
                       marker='s', edgecolors='black', linewidth=0.5)

            # Genome size bars (outer track)
            if size_track_enabled:
                size_val = row.get('size', np.nan)
                size_mbp = (size_val / 1e6) if pd.notna(size_val) else 0
                if max_size_mbp > 0 and size_mbp > 0:
                    size_length = (size_mbp / max_size_mbp) * size_track_height
                    if size_length > 0:
                        ax.bar(angle, size_length, width=size_bar_width,
                               bottom=size_base_radius, align='center',
                               color=self.colors['size_color'], edgecolor='black',
                               linewidth=0.3, alpha=0.85)

            # Optional species labels for small datasets only
            if show_tip_labels:
                label = row.get('label', row.get('bin_id', ''))
                formatted_label = format_species_name(label)
                rotation = np.degrees(angle)
                if 90 < rotation < 270:
                    rotation += 180
                    ha = 'right'
                else:
                    ha = 'left'

                ax.text(angle, label_radius, formatted_label,
                        rotation=rotation, ha=ha, va='center',
                        fontsize=7, alpha=0.9)

        if current_phylum is not None:
            phylum_segments.append((current_phylum, phylum_start_idx, n_bins))

        # Phylum color ring on the outermost track
        for phylum, start_idx, end_idx in phylum_segments:
            if end_idx <= start_idx:
                continue
            start_angle = angle_offset + start_idx * angle_step
            end_angle = angle_offset + end_idx * angle_step
            span = end_angle - start_angle
            mid_angle = start_angle + span / 2
            color = phylum_colors.get(phylum, '#cccccc')
            ax.bar(mid_angle, phylum_ring_width, width=span, bottom=phylum_ring_inner,
                   color=color, edgecolor='none', align='center', linewidth=0)

        # Customize
        ax.set_ylim(0, plot_max_radius)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.grid(False)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.spines['polar'].set_visible(False)
        
        # Legends
        if size_track_enabled:
            gap_center = 0.0  # gap is centred on 0 radians by construction
            self._draw_size_scale(ax, size_base_radius, size_track_height, max_size_mbp, label_theta=gap_center)

        self._add_circular_legend(fig, include_size_track=size_track_enabled)
        self._add_phylum_legend(fig, phylum_colors)
        
        plt.tight_layout(rect=[0.02, 0.08, 0.98, 0.98])
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close()

        return str(output_path)
    
    def _create_rectangular_tree(self, output_path, dpi=300):
        """Create rectangular tree for ≤30 bins."""
        n_bins = len(self.df)
        figsize = (12, max(8, n_bins * 0.3))

        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 4, width_ratios=[0.92, 0.10, 0.10, 0.52],
                              wspace=0.00, hspace=0)
        ax_tree = fig.add_subplot(gs[0])
        ax_comp = fig.add_subplot(gs[1])
        ax_cov  = fig.add_subplot(gs[2])
        ax_size = fig.add_subplot(gs[3])

        branch_end = 0.24
        label_x = 0.302
        y_spacing = 1.0

        tree_layout = None
        if self.tree:
            tree_layout = self._build_rectangular_tree_layout(branch_end, y_spacing=y_spacing)
            if tree_layout and 'bin_id' in self.df.columns:
                order_lookup = {bid: idx for idx, bid in enumerate(tree_layout['leaf_order'])}
                fallback = len(order_lookup)
                self.df['__tree_order'] = self.df['bin_id'].map(lambda x: order_lookup.get(x, fallback))
                self.df.sort_values('__tree_order', kind='mergesort', inplace=True)
                self.df.drop(columns='__tree_order', inplace=True)

        for col in ['phylum', 'class', 'order', 'family']:
            if col not in self.df.columns:
                self.df[col] = 'Unclassified' if col == 'phylum' else 'Unknown'

        if not tree_layout:
            if not self.bin_order:
                self.df = self.df.sort_values(['phylum', 'class', 'order', 'family'], na_position='last')

        phylum_colors = {}
        for i, phylum in enumerate(self.df['phylum'].unique()):
            phylum_colors[phylum] = self.colors['phylum_palette'][i % len(self.colors['phylum_palette'])]

        y_position = 0
        fallback_tree_y = 0
        phylum_y_positions = {}  # Track all y positions per phylum

        comp_entries = []
        cov_entries = []
        y_values = []

        if tree_layout:
            for edge in tree_layout['edges']:
                ax_tree.plot([edge['px'], edge['px']], [edge['py'], edge['cy']],
                             color='#8c8c8c', linewidth=1.2, alpha=0.8)
                ax_tree.plot([edge['px'], edge['cx']], [edge['cy'], edge['cy']],
                             color='#8c8c8c', linewidth=1.2, alpha=0.8)

        for _, row in self.df.iterrows():
            if tree_layout and 'bin_id' in row:
                y_val = tree_layout['leaf_y'].get(row['bin_id'])
                if y_val is None:
                    y_val = fallback_tree_y
                    fallback_tree_y += y_spacing
            else:
                y_val = y_position

            if not tree_layout:
                ax_tree.plot([0.045, branch_end], [y_val, y_val],
                             'k-', linewidth=1, alpha=0.5)

            phylum_value = row.get('phylum', 'Unknown')
            if phylum_value not in phylum_y_positions:
                phylum_y_positions[phylum_value] = []
            phylum_y_positions[phylum_value].append(y_val)

            if row.get('is_MAG', False):
                color = (self.colors['circular_mag_color']
                         if row.get('is_fully_circular', False)
                         else self.colors['mag_color'])
                ax_tree.scatter(branch_end, y_val, s=60, c=[color],
                                marker='o', edgecolors='none', zorder=6)

            label = row.get('label', row.get('bin_id', ''))
            formatted_label = format_species_name(label)
            ax_tree.text(label_x, y_val, formatted_label,
                         va='center', fontsize=10)

            completeness = row.get('Completeness', 0)
            comp_entries.append((y_val, completeness))

            coverage = row.get('mean_coverage', 1)
            if coverage > 0:
                log_cov = np.log10(coverage)
                cov_norm = min(max(log_cov / 3, 0), 1)
            else:
                cov_norm = 0
            cov_entries.append((y_val, cov_norm))
            size_mbp = row.get('size', 0) / 1e6
            ax_size.barh(
                y_val, size_mbp, height=0.6,
                color=self.colors['size_color'],
                edgecolor='black', linewidth=0.5, alpha=0.8
            )

            y_values.append(y_val)

            if not tree_layout:
                y_position += y_spacing

        # Draw phylum bars - handle non-contiguous groups
        # Group consecutive bins of the same phylum
        phylum_segments = []
        for phylum in phylum_y_positions.keys():
            positions = sorted(phylum_y_positions[phylum])
            if not positions:
                continue
            
            # Find contiguous segments
            segments = []
            current_segment = [positions[0]]
            
            for i in range(1, len(positions)):
                # Check if consecutive (within tolerance for y_spacing)
                if abs(positions[i] - positions[i-1]) <= y_spacing * 1.5:
                    current_segment.append(positions[i])
                else:
                    # Start new segment
                    segments.append(current_segment)
                    current_segment = [positions[i]]
            segments.append(current_segment)
            
            # Draw each segment
            for segment in segments:
                start_y = min(segment)
                end_y = max(segment)
                # For singletons, create a small bar centered on the bin
                if end_y <= start_y:
                    start_y = start_y - y_spacing * 0.4
                    end_y = end_y + y_spacing * 0.4
                phylum_segments.append((phylum, start_y, end_y))
        
        # Now draw all segments
        for phylum, start_y, end_y in phylum_segments:
            color = phylum_colors.get(phylum, '#cccccc')
            ax_tree.plot([0, 0], [start_y, end_y],
                         color=color, linewidth=4, alpha=0.8)
            label = str(phylum).replace('Candidatus', 'Ca.')
            mid_y = (start_y + end_y) / 2
            ax_tree.text(-0.025, mid_y, label,
                         ha='right', va='center', fontsize=9, fontweight='bold',
                         color=color)

        if tree_layout and not y_values:
            y_values = [0.0]

        y_min = min(y_values) - y_spacing if y_values else -1
        y_max = max(y_values) + y_spacing if y_values else (y_position if y_position else 1)

        ax_tree.set_xlim(-0.035, 0.69)
        ax_tree.set_ylim(y_min, y_max)
        ax_tree.axis('off')

        for ax, title in [
            (ax_comp, 'Comp\n(%)'),
            (ax_cov, 'Cov\n(log10)'),
        ]:
            ax.set_xlim(0, 1)
            ax.set_ylim(y_min, y_max)
            ax.axis('off')
            ax.set_title(title, fontsize=9, pad=10)

        fig.subplots_adjust(left=0.09, right=0.985, top=0.97, bottom=0.35, wspace=0.00)

        # Draw completeness/coverage squares using axis coordinates so they stay square
        y_range = y_max - y_min if y_max > y_min else 1
        row_height_frac = y_spacing / y_range
        square_height = row_height_frac * 0.78
        fig_w, fig_h = fig.get_size_inches()

        def _axis_height_width_ratio(axis):
            pos = axis.get_position()
            axis_height = pos.height * fig_h
            axis_width = pos.width * fig_w
            return axis_height / axis_width if axis_width else 1

        width_per_height = {
            ax_comp: _axis_height_width_ratio(ax_comp),
            ax_cov: _axis_height_width_ratio(ax_cov),
        }

        def _add_square(axis, x_center, y_val, color):
            y_frac = (y_val - y_min) / y_range
            axis_ratio = width_per_height[axis]
            square_width = square_height * axis_ratio
            axis.add_patch(Rectangle(
                (x_center - square_width / 2, y_frac - square_height / 2),
                square_width, square_height,
                transform=axis.transAxes,
                facecolor=color,
                edgecolor='black',
                linewidth=0.5,
                clip_on=False,
            ))

        for y_val, completeness in comp_entries:
            comp_color = self.colors['completeness_cmap'](completeness / 100)
            _add_square(ax_comp, 0.48, y_val, comp_color)

        for y_val, cov_norm in cov_entries:
            cov_color = self.colors['coverage_cmap'](cov_norm)
            _add_square(ax_cov, 0.52, y_val, cov_color)

        max_size = self.df['size'].max() / 1e6 if 'size' in self.df.columns else 10
        ax_size.set_xlim(0, max_size * 1.1)
        ax_size.set_ylim(y_min, y_max)
        ax_size.set_xlabel('Genome size (Mbp)', fontsize=10)
        #ax_size.set_title('Genome size', fontsize=10, pad=10)
        ax_size.set_yticks([])
        ax_size.spines['top'].set_visible(False)
        ax_size.spines['right'].set_visible(False)
        ax_size.spines['left'].set_visible(False)
        ax_size.grid(axis='x', alpha=0.3, linestyle='--')


        self._add_rectangular_legend(fig)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close()

        return str(output_path)

    
    def _add_circular_legend(self, fig, include_size_track=False):
        """Add legend for circular layout."""
        legend_elements = []
        
        # MAG indicators
        from matplotlib.lines import Line2D
        legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                     markerfacecolor=self.colors['mag_color'],
                                     markersize=8, label='MAG'))
        legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                     markerfacecolor=self.colors['circular_mag_color'],
                                     markersize=8, label='Circular MAG'))
        
        # Quality indicators
        legend_elements.append(mpatches.Patch(facecolor=self.colors['completeness_cmap'](1.0),
                                             label='High completeness'))
        legend_elements.append(mpatches.Patch(facecolor=self.colors['coverage_cmap'](1.0),
                                             label='High coverage'))

        if include_size_track:
            legend_elements.append(mpatches.Patch(facecolor=self.colors['size_color'],
                                                 label='Genome size (Mbp)'))
        
        plt.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True)


    def _add_phylum_legend(self, fig, phylum_colors):
        """Display phylum colors separately from the outer ring."""
        if not phylum_colors:
            return

        legend_handles = []
        for phylum, color in phylum_colors.items():
            if pd.notna(phylum) and str(phylum).strip():
                label = str(phylum)
            else:
                label = 'Unclassified'
            label = label.replace('Candidatus', 'Ca.')
            legend_handles.append(mpatches.Patch(facecolor=color, label=label))

        if not legend_handles:
            return

        n_cols = min(5, len(legend_handles))
        fig.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, 0.02),
                   ncol=n_cols, frameon=False, columnspacing=0.8, handletextpad=0.4,
                   fontsize=10)


    def _draw_size_scale(self, ax, base_radius, track_height, max_size_mbp, label_theta=0.0):
        """Overlay dashed rings and labels to provide a genome-size scale."""
        if max_size_mbp <= 0:
            return

        theta = np.linspace(0, 2 * np.pi, 720)
        fractions = [0.25, 0.5, 0.75, 1.0]

        for frac in fractions:
            radius = base_radius + track_height * frac
            ax.plot(theta, np.full_like(theta, radius), color='#bdbdbd', linestyle='--',
                    linewidth=0.4, alpha=0.7, zorder=1)
            label = f"{max_size_mbp * frac:.1f} Mbp"
            ax.text(label_theta, radius, label, fontsize=8, color='#555555',
                    ha='right', va='center')


    def _add_rectangular_legend(self, fig):
        """Add legend for rectangular layout below the axes."""
        legend_base = 0.02
        footer_top = 0.40  # matches subplots_adjust bottom
        footer_height = footer_top - legend_base

        # MAG indicators block (upper part of footer)
        mag_block_y = legend_base + 0.77 * footer_height
        fig.text(0.10, mag_block_y + 0.02, 'MAG indicators:',
                 fontsize=10, fontweight='bold', va='bottom')
        mag_ax = fig.add_axes([0.14, mag_block_y - 0.02, 0.15, 0.045])
        mag_ax.axis('off')
        mag_ax.scatter(0.15, 0.5, s=60, c=[self.colors['mag_color']],
                       marker='o', edgecolors='black', linewidth=0.5)
        mag_ax.text(0.25, 0.5, 'MAG', fontsize=9, va='center')
        mag_ax.scatter(0.65, 0.5, s=60, c=[self.colors['circular_mag_color']],
                       marker='o', edgecolors='black', linewidth=0.5)
        mag_ax.text(0.75, 0.5, 'Circular MAG', fontsize=9, va='center')
        mag_ax.set_xlim(0, 1)
        mag_ax.set_ylim(0, 1)

        # Completeness gradient stacked beneath MAG block
        comp_block_y = legend_base + 0.60 * footer_height
        fig.text(0.10, comp_block_y + 0.025, 'Completeness (%):', fontsize=10,
                 fontweight='bold', va='bottom')
        comp_ax = fig.add_axes([0.14, comp_block_y - 0.01, 0.15, 0.02])
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        comp_ax.imshow(gradient, aspect='auto', cmap=self.colors['completeness_cmap'])
        comp_ax.set_xticks([0, 128, 256])
        comp_ax.set_xticklabels(['0%', '50%', '100%'], fontsize=8)
        comp_ax.set_yticks([])
        for spine in comp_ax.spines.values():
            spine.set_visible(False)

        # Coverage gradient stacked at the bottom
        cov_block_y = legend_base + 0.35 * footer_height
        fig.text(0.10, cov_block_y + 0.025, 'Coverage (log10):', fontsize=10,
                 fontweight='bold', va='bottom')
        cov_ax = fig.add_axes([0.14, cov_block_y - 0.01, 0.15, 0.02])
        cov_ax.set_yticks([])
        grad = np.linspace(0, 1, 256).reshape(1, -1)
        cov_ax.imshow(grad, aspect='auto', cmap=self.colors['coverage_cmap'])
        max_cov = max(float(pd.to_numeric(self.df['mean_coverage'], errors='coerce').max()), 1.0) if 'mean_coverage' in self.df.columns else 1.0
        max_log = np.ceil(np.log10(max_cov))
        cov_ax.set_xticks([0, 128, 256])
        cov_ax.set_xticklabels([r'$10^{0}$', r'$10^{%d}$' % int(max_log//2), 
                                r'$10^{%d}$' % int(max_log)], fontsize=8)
        cov_ax.set_yticks([])
        for spine in cov_ax.spines.values():
            spine.set_visible(False)


def main():

    parser = argparse.ArgumentParser(description="Generate metagenome bin tree figure")
    parser.add_argument("--csv", required=True, help="Input bin_data.csv file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--dpi", type=int, default=300, help="Figure resolution (DPI)")
    parser.add_argument("--filename", default="metagenome_tree.png", help="Output filename")
    parser.add_argument("--order-json", help="Optional taxonomy_tree.json to reuse bin ordering")
    parser.add_argument("--tree-newick", help="Explicit Newick file to use for branch structure (overrides JSON entry)")
    
    args = parser.parse_args()
    
    # Setup paths
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / args.filename
    
    # Generate figure
    bin_order = None
    tree_path = None
    if args.order_json:
        bin_order, tree_path = load_order_metadata(args.order_json)
    if args.tree_newick:
        tree_path = Path(args.tree_newick)

    tree_fig = MetagenomeTreeFigure(args.csv, bin_order=bin_order, tree_path=tree_path)
    n_bins = tree_fig.load_and_process()
    
    if n_bins > 30:
        print(f"Creating circular tree for {n_bins} bins...")
    else:
        print(f"Creating rectangular tree for {n_bins} bins...")
    
    tree_fig.create_figure(output_path, dpi=args.dpi)
    print(f"Tree figure saved to: {output_path}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
