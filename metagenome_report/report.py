"""Generate metagenome summary table and context dict from bin_data.csv."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

try:
    from num2words import num2words
    def _text_num(n: int) -> str:
        return num2words(n) if n <= 10 else str(n)
except ImportError:
    def _text_num(n: int) -> str:
        return str(n)


# ─── column resolution ────────────────────────────────────────────────────────

class _Cols:
    def __init__(self, df: pd.DataFrame):
        self.map = {c.lower(): c for c in df.columns}

    def get(self, *cands: str) -> Optional[str]:
        for c in cands:
            if c in self.map:
                return self.map[c]
        return None


def _resolve_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    c = _Cols(df)
    return {
        "assembler":        c.get("assembler"),
        "binner":           c.get("binning_program", "binner", "binning"),
        "refiner":          c.get("refining_program", "bin_refinement", "refiner"),
        "bin_id":           c.get("bin_id", "bin", "bin_name"),
        "size":             c.get("size", "length", "span", "length_bp", "total_length"),
        "contigs":          c.get("contigs", "n_contigs", "num_contigs"),
        "circular":         c.get("circular", "is_circular", "circularised", "circularized"),
        "bin_type":         c.get("bin_type", "bin classification", "bintype"),
        "trna_unique":      c.get("unique_trnas", "n_unique_trna", "trna_unique"),
        "rrna_5s":          c.get("rrna_5s", "5s_rrna", "has_5s"),
        "rrna_16s":         c.get("rrna_16s", "16s_rrna", "has_16s"),
        "rrna_23s":         c.get("rrna_23s", "23s_rrna", "has_23s"),
        "completeness":     c.get("completeness", "checkm_completeness", "Completeness"),
        "contamination":    c.get("contamination", "checkm_contamination", "Contamination"),
        "taxonomy":         c.get("classification", "gtdb_taxonomy", "gtdbtk_taxonomy"),
        "ncbi_taxonomy":    c.get("ncbi_classification"),
        "drep":             c.get("drep", "drep_cluster"),
        "checkm_version":   c.get("checkm_version"),
        "checkm_db":        c.get("checkm_db", "checkm_database"),
        "gtdbtk_version":   c.get("gtdbtk_version"),
        "gtdb_release":     c.get("gtdb_release", "gtdb_version"),
        "mean_coverage":    c.get("mean_coverage", "coverage", "avg_coverage"),
    }


# ─── normalisation helpers ────────────────────────────────────────────────────

_KNOWN_BINNERS: List[str] = [
    "metabat2", "maxbin2", "concoct", "vamb", "semibin", "rosella",
    "bin3c", "graphbin", "binsanity", "autometa",
]


def _norm_assembler(s: Optional[str]) -> Optional[str]:
    if not s:
        return None
    v = s.strip().lower()
    if v in {"meta-mdbg", "metamdbg", "meta_mdbg"}:
        return "metamdbg"
    return v


def _norm_binner(s: Optional[str]) -> Optional[str]:
    if not s:
        return None
    v = s.strip().lower()
    for k in _KNOWN_BINNERS:
        if k in v:
            return k
    if "magscot" in v or "dastool" in v:
        return None
    toks = re.split(r"[^a-z0-9\+._-]+", v)
    return toks[0] if toks and toks[0] else None


def _norm_refiner(s: Optional[str]) -> Optional[str]:
    if not s:
        return None
    v = s.strip().lower()
    if "magscot" in v:
        return "magscot"
    if "dastool" in v or "das tool" in v:
        return "dastool"
    return v


def _is_true(v: Any) -> bool:
    if v is None:
        return False
    if isinstance(v, str):
        return v.strip().lower() in {"1", "true", "t", "yes", "y", "circular"}
    try:
        return float(v) >= 1.0
    except Exception:
        return False


def _as_float(v: Any) -> float:
    try:
        return float(v)
    except Exception:
        return float("nan")


def _extract_phylum(tax: Any) -> Optional[str]:
    if tax is None or (isinstance(tax, float) and np.isnan(tax)):
        return None
    m = re.search(r"p__([A-Za-z0-9_\-]+)", str(tax))
    return m.group(1) if m else None


def _extract_domain(tax: Any) -> Optional[str]:
    if tax is None or (isinstance(tax, float) and np.isnan(tax)):
        return None
    m = re.search(r"d__([A-Za-z0-9_]+)", str(tax))
    return m.group(1) if m else None


# ─── MAG classification ────────────────────────────────────────────────────────

def _has_rrna_operon(row: pd.Series, cols: Dict[str, Optional[str]]) -> Optional[bool]:
    r5  = str(row[cols["rrna_5s"]]).upper() == "Y"  if cols["rrna_5s"]  else None
    r16 = str(row[cols["rrna_16s"]]).upper() == "Y" if cols["rrna_16s"] else None
    r23 = str(row[cols["rrna_23s"]]).upper() == "Y" if cols["rrna_23s"] else None
    if None in (r5, r16, r23):
        return None
    return bool(r5 and r16 and r23)


def _is_mag(row: pd.Series, cols: Dict[str, Optional[str]]) -> bool:
    if cols.get("bin_type"):
        bt = str(row[cols["bin_type"]]).strip().lower()
        if bt:
            return "mag" in bt

    comp = _as_float(row[cols["completeness"]]) if cols["completeness"] else float("nan")
    cont = _as_float(row[cols["contamination"]]) if cols["contamination"] else float("nan")

    if np.isnan(comp) or np.isnan(cont) or cont > 5:
        return False

    rrna_ok = _has_rrna_operon(row, cols)
    if rrna_ok is not None and not rrna_ok:
        return False

    trna_ok: Optional[bool] = None
    if cols["trna_unique"]:
        try:
            trna_ok = int(row[cols["trna_unique"]]) >= 18
        except Exception:
            trna_ok = None
    if trna_ok is not None and not trna_ok:
        return False

    fully_circular = False
    if cols["contigs"] and cols["circular"]:
        try:
            nc = float(row[cols["contigs"]])
            ncirc = float(row[cols["circular"]])
            fully_circular = (nc == 1 and ncirc == 1) or (nc > 1 and nc == ncirc)
        except Exception:
            pass
    elif cols["circular"]:
        fully_circular = _is_true(row[cols["circular"]])

    return (comp >= 90.0) or (comp >= 50.0 and fully_circular)


# ─── public API ───────────────────────────────────────────────────────────────

def build_context(csv_path: Path) -> Dict[str, Any]:
    """Return a metagenome context dict for a given bin_data CSV."""
    df = pd.read_csv(csv_path)
    cols = _resolve_columns(df)

    def _first(colname: Optional[str]) -> Optional[str]:
        return df[colname].dropna().astype(str).iloc[0] if colname and df[colname].notna().any() else None

    assembler = _norm_assembler(_first(cols["assembler"]))
    binners_raw = df[cols["binner"]].dropna().astype(str).unique().tolist() if cols["binner"] else []
    binners = sorted({b for b in (_norm_binner(x) for x in binners_raw) if b})
    refiner = _norm_refiner(_first(cols["refiner"]))

    tax_col = cols["taxonomy"] or cols["ncbi_taxonomy"]
    if tax_col:
        df["domain"] = df[tax_col].apply(_extract_domain)
        df["phylum"] = df[tax_col].apply(_extract_phylum)
        has_archaea = (df["domain"] == "Archaea").any()
        num_phyla = int(df["phylum"].nunique(dropna=True))
    else:
        has_archaea = False
        num_phyla = 0

    df["is_MAG"] = df.apply(lambda r: _is_mag(r, cols), axis=1)

    total_bins = len(df)

    # circular MAG count
    if cols["contigs"] and cols["circular"]:
        nc = pd.to_numeric(df[cols["contigs"]], errors="coerce")
        ncirc = pd.to_numeric(df[cols["circular"]], errors="coerce")
        fully_circular = (((nc == 1) & (ncirc == 1)) | ((nc > 1) & (nc == ncirc))).fillna(False)
        num_circular_mags = int((df["is_MAG"] & fully_circular).sum())
    elif cols["circular"]:
        circ_series = df[cols["circular"]].astype(str).str.lower().isin(
            ["1", "y", "yes", "true", "t", "circular"]
        )
        num_circular_mags = int((df["is_MAG"] & circ_series).sum())
    else:
        num_circular_mags = 0

    num_mags = int(df["is_MAG"].sum())

    if cols["size"] and df[cols["size"]].notna().any():
        sizes_mbp = pd.to_numeric(df[cols["size"]], errors="coerce") / 1e6
        min_size_mbp  = float(np.nanmin(sizes_mbp))
        max_size_mbp  = float(np.nanmax(sizes_mbp))
        mean_size_mbp = float(np.nanmean(sizes_mbp))
        std_size_mbp  = float(np.nanstd(sizes_mbp))
    else:
        min_size_mbp = max_size_mbp = mean_size_mbp = std_size_mbp = float("nan")

    if cols["completeness"] and cols["contamination"]:
        comp = pd.to_numeric(df[cols["completeness"]], errors="coerce")
        cont = pd.to_numeric(df[cols["contamination"]], errors="coerce")
        mean_completeness  = float(np.nanmean(comp))
        std_completeness   = float(np.nanstd(comp))
        mean_contamination = float(np.nanmean(cont))
        std_contamination  = float(np.nanstd(cont))
    else:
        mean_completeness = std_completeness = mean_contamination = std_contamination = float("nan")

    def _r(v: float) -> Optional[float]:
        return round(v, 6) if np.isfinite(v) else None

    return {
        "assembler":           assembler,
        "binners":             binners,
        "multiple_binners":    len(binners) >= 2,
        "refiner":             refiner,
        "checkm_version":      _first(cols["checkm_version"]),
        "checkm_db":           _first(cols["checkm_db"]),
        "gtdbtk_version":      _first(cols["gtdbtk_version"]),
        "gtdb_release":        _first(cols["gtdb_release"]),
        "use_drep":            bool(cols["drep"]),
        "drep_threshold":      None,
        "filter_domains":      False,
        "mag_contamination":   5,
        "min_trnas":           18,
        "high_completeness":   90,
        "medium_completeness": 50,
        "min_bin_completeness": 50,
        "max_bin_contamination": 10,
        "total_bins":          total_bins,
        "num_mags":            num_mags,
        "num_circular_mags":   num_circular_mags,
        "num_phyla":           num_phyla,
        "total_bins_text":     _text_num(total_bins),
        "num_mags_text":       _text_num(num_mags),
        "num_circular_mags_text": _text_num(num_circular_mags),
        "num_phyla_text":      _text_num(num_phyla),
        "has_archaea":         bool(has_archaea),
        "min_size_mbp":        _r(min_size_mbp),
        "max_size_mbp":        _r(max_size_mbp),
        "mean_size_mbp":       _r(mean_size_mbp),
        "std_size_mbp":        _r(std_size_mbp),
        "mean_completeness":   _r(mean_completeness),
        "std_completeness":    _r(std_completeness),
        "mean_contamination":  _r(mean_contamination),
        "std_contamination":   _r(std_contamination),
        "visualization_method": "custom",
        "has_figure":          False,
        "figure_relpath":      None,
    }


def write_table(csv_path: Path, outdir: Path) -> tuple[Path, Path]:
    """Write markdown and CSV summary tables of bins. Returns (md_path, csv_path)."""
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_path)
    cols = _resolve_columns(df)

    df["is_MAG"] = df.apply(lambda r: _is_mag(r, cols), axis=1)

    tax_col = cols["taxonomy"] or cols["ncbi_taxonomy"]
    if tax_col:
        df["phylum"] = df[tax_col].apply(_extract_phylum)

    keep_candidates = [
        cols["bin_id"], "phylum", tax_col,
        cols["completeness"], cols["contamination"],
        cols["circular"], cols["trna_unique"],
        cols["rrna_5s"], cols["rrna_16s"], cols["rrna_23s"],
        "is_MAG",
    ]
    keep = [c for c in keep_candidates if c and c in df.columns]
    # deduplicate preserving order
    seen: set = set()
    keep = [c for c in keep if not (c in seen or seen.add(c))]  # type: ignore[func-returns-value]

    summary = df[keep]
    md_path = outdir / "metagenome_bins_table.md"
    md_path.write_text(summary.to_markdown(index=False), encoding="utf-8")

    csv_out = outdir / "metagenome_bins_table.csv"
    summary.to_csv(csv_out, index=False)

    return md_path, csv_out
