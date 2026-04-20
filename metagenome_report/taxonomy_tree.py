"""Build a rank-aware ete3 taxonomy tree from bin_data.csv."""

from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

try:
    from ete3 import NCBITaxa, Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    Tree = None
    NCBITaxa = None


RANK_ORDER: Tuple[str, ...] = (
    "domain", "phylum", "class", "order", "family", "genus", "species",
)

PREFIX_TO_RANK = {
    "d__": "domain",
    "k__": "domain",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}

RANK_TO_PREFIX = {v: k for k, v in PREFIX_TO_RANK.items() if len(k) == 3}


def parse_classification(classification: str) -> Dict[str, str]:
    if not classification:
        return {}
    ranks: Dict[str, str] = {}
    for token in str(classification).split(";"):
        token = token.strip()
        if not token or "__" not in token:
            continue
        prefix, _, remainder = token.partition("__")
        prefix = f"{prefix}__"
        rank = PREFIX_TO_RANK.get(prefix)
        if not rank:
            continue
        name = remainder.strip()
        if not name or name == "Unknown":
            continue
        ranks[rank] = name
    return ranks


def _sanitize_label(value: str) -> str:
    return (
        value.replace(" ", "_").replace("/", "_")
        .replace("(", "").replace(")", "").replace(",", "")
    )


class TaxonomyTreeBuilder:
    def __init__(
        self,
        dataframe: pd.DataFrame,
        taxonomy_columns: Iterable[str] = ("ncbi_classification", "classification"),
        use_taxid_lookup: bool = False,
        ncbi_db: Optional[str] = None,
        rank_order: Tuple[str, ...] = RANK_ORDER,
    ) -> None:
        if not ETE3_AVAILABLE:
            raise ImportError("ete3 is required. Install with: pip install ete3")
        self.df = dataframe.copy()
        self.taxonomy_columns = tuple(taxonomy_columns)
        self.rank_order = rank_order
        self.ncbi: Optional[object] = None
        if use_taxid_lookup:
            self.ncbi = NCBITaxa(dbfile=ncbi_db) if ncbi_db else NCBITaxa()
        self.node_cache: dict = {}

    def build(self) -> Tuple[object, List[str]]:
        root = Tree(name="root", dist=0.0)
        self.node_cache = {tuple(): root}

        for _, row in self.df.iterrows():
            bin_id = str(row.get("bin_id", "")).strip()
            if not bin_id:
                continue
            taxonomy = self._extract_taxonomy(row)
            node = root
            path: List[Tuple[str, str]] = []

            for rank in self.rank_order:
                label = taxonomy.get(rank)
                if not label:
                    continue
                safe = f"{RANK_TO_PREFIX.get(rank, rank[0] + '__')}{_sanitize_label(label)}"
                path.append((rank, safe))
                cache_key = tuple(path)
                if cache_key not in self.node_cache:
                    child = node.add_child(name=safe, dist=1.0)
                    child.add_features(rank=rank, label=label)
                    self.node_cache[cache_key] = child
                node = self.node_cache[cache_key]

            leaf_key = tuple(path + [("bin", bin_id)])
            if leaf_key in self.node_cache:
                continue
            leaf = node.add_child(name=bin_id, dist=1.0)
            leaf.add_features(rank="bin", label=bin_id)
            self.node_cache[leaf_key] = leaf

        root.sort_descendants()
        leaf_order = [leaf.name for leaf in root.iter_leaves()]
        return root, leaf_order

    def _extract_taxonomy(self, row: pd.Series) -> Dict[str, str]:
        for column in self.taxonomy_columns:
            if column in row and pd.notna(row[column]):
                parsed = parse_classification(row[column])
                if parsed:
                    return parsed
        if self.ncbi is not None and "taxon_id" in row and pd.notna(row["taxon_id"]):
            try:
                taxid = int(row["taxon_id"])
            except (TypeError, ValueError):
                taxid = None
            if taxid:
                return self._taxonomy_from_taxid(taxid)
        return {}

    def _taxonomy_from_taxid(self, taxid: int) -> Dict[str, str]:
        try:
            lineage = self.ncbi.get_lineage(taxid)
        except Exception:
            return {}
        names = self.ncbi.get_taxid_translator(lineage)
        ranks = self.ncbi.get_rank(lineage)
        mapping: Dict[str, str] = {}
        for tid in lineage:
            rank = ranks.get(tid)
            # NCBI uses "superkingdom" where we use "domain"
            if rank == "superkingdom":
                rank = "domain"
            if rank in self.rank_order:
                mapping[rank] = names.get(tid, str(tid))
        return mapping


def build_taxonomy_tree(
    csv_path: Path,
    outdir: Path,
    taxonomy_columns: Iterable[str] = ("ncbi_classification", "classification"),
    use_taxid_lookup: bool = False,
    ncbi_db: Optional[str] = None,
    newick_name: str = "taxonomy_tree.nwk",
    metadata_name: str = "taxonomy_tree.json",
) -> Tuple[Path, Path]:
    """Build taxonomy tree from CSV; write Newick + JSON metadata. Returns (nwk_path, json_path)."""
    df = pd.read_csv(csv_path)
    if "bin_id" not in df.columns:
        raise ValueError("CSV must contain a 'bin_id' column")
    df = df.drop_duplicates(subset="bin_id", keep="first")

    outdir.mkdir(parents=True, exist_ok=True)
    builder = TaxonomyTreeBuilder(
        df,
        taxonomy_columns=list(taxonomy_columns),
        use_taxid_lookup=use_taxid_lookup,
        ncbi_db=ncbi_db,
    )
    tree, leaf_order = builder.build()

    newick_path = outdir / newick_name
    tree.write(outfile=str(newick_path), format=1)

    meta = OrderedDict(
        tree_newick=str(newick_path),
        bin_order=leaf_order,
        rank_order=list(builder.rank_order),
        taxonomy_columns=list(taxonomy_columns),
        used_taxid_lookup=bool(builder.ncbi),
    )
    metadata_path = outdir / metadata_name
    metadata_path.write_text(json.dumps(meta, indent=2))
    return newick_path, metadata_path
