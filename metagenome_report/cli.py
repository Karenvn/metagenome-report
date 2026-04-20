#!/usr/bin/env python3
"""
metagenome-report: generate bin table and tree figure from gn_assets data.

Usage examples
--------------
# By tolid (auto-finds CSV in ~/gn_assets/metagenomes/<tolid>/bin_data.csv)
metagenome-report --tolid glLicPygm2

# Explicit CSV + explicit output dir
metagenome-report --csv /path/to/bin_data.csv --outdir /path/to/output

# Also build taxonomy tree first (requires ete3)
metagenome-report --tolid glLicPygm2 --build-tree

# Use pre-built taxonomy tree JSON
metagenome-report --tolid glLicPygm2 --order-json /path/to/taxonomy_tree.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

GN_ASSETS = Path.home() / "gn_assets"
METAGENOMES_DIR = GN_ASSETS / "metagenomes"
FIGS_DIR = GN_ASSETS / "metagenome_figs"


def _find_csv(tolid: str) -> Path:
    candidate = METAGENOMES_DIR / tolid / "bin_data.csv"
    if not candidate.is_file():
        raise FileNotFoundError(
            f"No bin_data.csv found for tolid '{tolid}' at {candidate}"
        )
    return candidate


def _default_outdir(tolid: Optional[str], csv_path: Path) -> Path:
    if tolid:
        return FIGS_DIR / tolid
    return FIGS_DIR / csv_path.parent.name


def run(
    csv_path: Path,
    outdir: Path,
    tolid: Optional[str] = None,
    build_tree: bool = False,
    order_json: Optional[Path] = None,
    tree_newick: Optional[Path] = None,
    use_taxid: bool = False,
    ncbi_db: Optional[str] = None,
    dpi: int = 300,
    figure_filename: str = "metagenome_tree.png",
    skip_figure: bool = False,
) -> None:
    from metagenome_report.report import build_context, write_table
    from metagenome_report.tree_figure import MetagenomeTreeFigure, load_order_metadata

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {outdir}")

    # 1. Optionally build taxonomy tree
    tree_assets_dir = outdir / "tree_assets"
    auto_json = tree_assets_dir / "taxonomy_tree.json"

    if build_tree:
        try:
            from metagenome_report.taxonomy_tree import build_taxonomy_tree
            print("Building taxonomy tree …")
            nwk, meta_json = build_taxonomy_tree(
                csv_path=csv_path,
                outdir=tree_assets_dir,
                use_taxid_lookup=use_taxid,
                ncbi_db=ncbi_db,
            )
            print(f"  Newick:   {nwk}")
            print(f"  Metadata: {meta_json}")
            order_json = meta_json
        except ImportError as exc:
            print(f"[WARN] ete3 not available, skipping tree build: {exc}", file=sys.stderr)
    elif order_json is None and auto_json.is_file():
        order_json = auto_json
        print(f"Found existing taxonomy metadata: {order_json}")

    # 2. Generate figure
    if not skip_figure:
        bin_order = None
        resolved_tree: Optional[Path] = tree_newick
        if order_json and Path(order_json).is_file():
            try:
                bin_order, resolved_tree = load_order_metadata(str(order_json))
                print(f"Using taxonomy order from {order_json}")
            except Exception as exc:
                print(f"[WARN] Could not read order JSON: {exc}", file=sys.stderr)

        if tree_newick:
            resolved_tree = tree_newick

        figure_path = outdir / figure_filename
        print(f"Generating tree figure ({figure_path.name}) …")
        try:
            fig = MetagenomeTreeFigure(
                str(csv_path),
                bin_order=bin_order,
                tree_path=resolved_tree,
            )
            n = fig.load_and_process()
            fig.create_figure(figure_path, dpi=dpi)
            print(f"  Figure saved: {figure_path}  ({n} bins)")
        except Exception as exc:
            print(f"[ERROR] Figure generation failed: {exc}", file=sys.stderr)

    # 3. Write markdown table
    table_path = write_table(csv_path, outdir)
    print(f"  Table saved:  {table_path}")

    # 4. Write context JSON
    ctx = build_context(csv_path)
    ctx["has_figure"] = not skip_figure and (outdir / figure_filename).is_file()
    ctx["figure_relpath"] = figure_filename if ctx["has_figure"] else None
    if tolid:
        ctx["tolid"] = tolid
    ctx_path = outdir / "metagenome_context.json"
    ctx_path.write_text(json.dumps(ctx, indent=2, default=str), encoding="utf-8")
    print(f"  Context JSON: {ctx_path}")


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate metagenome bin table and tree figure from gn_assets data."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--tolid", help="ToLID (auto-finds CSV in ~/gn_assets/metagenomes/)")
    src.add_argument("--csv", type=Path, help="Explicit path to bin_data.csv")

    p.add_argument("--outdir", type=Path,
                   help="Output directory (default: ~/gn_assets/metagenome_figs/<tolid>)")
    p.add_argument("--build-tree", action="store_true",
                   help="Build taxonomy tree with ete3 before plotting (writes tree_assets/)")
    p.add_argument("--order-json", type=Path,
                   help="taxonomy_tree.json from a prior --build-tree run")
    p.add_argument("--tree-newick", type=Path,
                   help="Explicit Newick file to use for branch structure")
    p.add_argument("--use-taxid", action="store_true",
                   help="Allow NCBI taxon_id lookups via ete3 (used with --build-tree)")
    p.add_argument("--ncbi-db", help="Path to ete3 NCBI taxonomy SQLite database")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300)")
    p.add_argument("--filename", default="metagenome_tree.png",
                   help="Figure filename (default: metagenome_tree.png)")
    p.add_argument("--skip-figure", action="store_true",
                   help="Only write table and context JSON, skip figure generation")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    tolid: Optional[str] = args.tolid
    if tolid:
        csv_path = _find_csv(tolid)
    else:
        csv_path = args.csv.expanduser().resolve()
        if not csv_path.is_file():
            print(f"[ERROR] CSV not found: {csv_path}", file=sys.stderr)
            return 1

    outdir = args.outdir or _default_outdir(tolid, csv_path)
    outdir = outdir.expanduser().resolve()

    run(
        csv_path=csv_path,
        outdir=outdir,
        tolid=tolid,
        build_tree=args.build_tree,
        order_json=args.order_json,
        tree_newick=args.tree_newick,
        use_taxid=args.use_taxid,
        ncbi_db=args.ncbi_db,
        dpi=args.dpi,
        figure_filename=args.filename,
        skip_figure=args.skip_figure,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
