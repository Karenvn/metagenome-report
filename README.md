# metagenome-report

Generate a metagenome bin summary table and annotated tree figure from
`bin_data.csv` files produced by the Tree of Life assembly/binning/QC
pipeline.  Outputs land in `~/gn_assets/metagenome_figs/<tolid>/` for
direct use in genome notes.


## Requirements

```
pandas  numpy  matplotlib  tabulate  num2words
```

Optional (for taxonomy-based tree layout):

```
ete3
```

Install the package and its core dependencies:

```bash
pip install -e .
# with ete3 support:
pip install -e ".[tree]"
```


## Input

Each sample is expected to have a `bin_data.csv` under:

```
~/gn_assets/metagenomes/<tolid>/bin_data.csv
```

The CSV must contain a `bin_id` column.  The tool recognises flexible
column naming (e.g. `Completeness` or `completeness`, `classification`
or `gtdb_taxonomy`, etc.) so files from different pipeline versions
should work without modification.


## Output

All outputs are written to `~/gn_assets/metagenome_figs/<tolid>/`
(override with `--outdir`):

| File | Description |
|---|---|
| `metagenome_tree.png` | Annotated bin tree figure (rectangular ≤50 bins, circular >50 bins) |
| `metagenome_bins_table.md` | Markdown summary table of bins |
| `metagenome_bins_table.csv` | Same table as CSV |
| `metagenome_context.json` | Context dict for genome note Jinja2 templates |
| `tree_assets/taxonomy_tree.nwk` | Newick tree (written when `--build-tree` is used) |
| `tree_assets/taxonomy_tree.json` | Bin order + rank metadata for the Newick tree |


## Usage

```bash
# Minimal — auto-finds CSV from tolid
metagenome-report --tolid glLicPygm2

# Explicit CSV path
metagenome-report --csv /path/to/bin_data.csv --outdir /path/to/output

# Build a GTDB-based taxonomy tree first, then plot with that topology
metagenome-report --tolid glLicPygm2 --build-tree

# Re-use a previously built taxonomy tree
metagenome-report --tolid glLicPygm2 --order-json ~/gn_assets/metagenome_figs/glLicPygm2/tree_assets/taxonomy_tree.json

# Allow NCBI taxon_id lookups when building the tree (needs ete3 + NCBI db)
metagenome-report --tolid glLicPygm2 --build-tree --use-taxid --ncbi-db ~/ncbi_taxadb.sqlite

# Only write table + context JSON, skip figure
metagenome-report --tolid glLicPygm2 --skip-figure

# Custom DPI and filename
metagenome-report --tolid glLicPygm2 --dpi 600 --filename metagenome_tree_hires.png
```

Run as a module if the entry point is not on your PATH:

```bash
python -m metagenome_report.cli --tolid glLicPygm2
```


## Plot appearance

**Rectangular layout** (≤50 bins): tree on the left, with coloured
phylum bars; completeness (blue gradient) and coverage (red gradient)
squares; genome size bars on the right.

**Circular layout** (>50 bins): bins arranged around a polar axis,
with an outer phylum colour ring, completeness/coverage scatter tracks,
and genome size bars.  Tip labels are shown only for ≤30 bins.

MAGs are marked with a filled circle: grey for standard MAGs, black for
fully circular MAGs (all contigs circular).

The figure uses Open Sans if available (`GENOMENOTES_FONT` env var can
point to an explicit `.ttf` path).


## Two-step workflow (taxonomy-ordered tree)

For a tree whose topology reflects GTDB classification, run
`--build-tree` once to produce the Newick file, then reuse it:

```bash
# Step 1 — build tree (only needed once per sample)
metagenome-report --tolid odAioCras1 --build-tree --use-taxid --ncbi-db ~/ncbi_taxadb.sqlite

# Step 2 — regenerate figure at higher DPI without rebuilding
metagenome-report --tolid odAioCras1 --dpi 600
```

The second run picks up `tree_assets/taxonomy_tree.json` automatically
if it exists.


## Batch workflow

Pass multiple tolids directly:

```bash
metagenome-report --tolids glLicPygm2 glLicPygm3 odAioCras1
```

Or supply a text file with one tolid per line (lines starting with `#` are ignored):

```bash
metagenome-report --tolids-file my_tolids.txt
```

Example `my_tolids.txt`:
```
# sponges
odAioCras1
odCarFoli1
# lichens
glLicPygm2
glLicPygm3
```

All other flags (`--build-tree`, `--dpi`, `--skip-figure`, etc.) apply to every sample in the batch. Each sample writes to its own `~/gn_assets/metagenome_figs/<tolid>/` directory. Samples that fail are reported at the end and do not stop the rest of the batch.


## Files

```
metagenome_report/
  cli.py            Command-line entry point
  report.py         build_context() and write_table()
  tree_figure.py    MetagenomeTreeFigure (rectangular + circular layouts)
  taxonomy_tree.py  TaxonomyTreeBuilder using ete3
pyproject.toml
```
