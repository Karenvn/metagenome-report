# merian-busco-painter

Utilities for plotting Merian element assignments from BUSCO `full_table.tsv` files,
adapted for genome notes workflows.

Merian elements are 32 ancestral chromosome building blocks found in the genomes of
moths and butterflies (Lepidoptera), which have remained largely intact for over
250 million years (Wright *et al*. 2024, doi:[10.1038/s41559-024-02329-4](https://doi.org/10.1038/s41559-024-02329-4)).

This code is derived from
[`charlottewright/lep_busco_painter`](https://github.com/charlottewright/lep_busco_painter)
and uses the same Merian reference table. The main changes in this version are:

- a Python plotting script (`plot_buscopainter.py`) instead of the original R plotter
- inclusion of duplicated BUSCO hits in the location table and plot
- chromosome length lookup from the NCBI Datasets API `sequence_reports` endpoint
- modifications to the appearance of the Merian plot for use in genome note publications
- automatic multi-column plotting for assemblies with more than 40 chromosomes
- a batch wrapper for running multiple ToLIDs from genome notes working directories.

## Important assumption

The bundled `Merian_elements_full_table.tsv` reference table is based on BUSCO
`odb10`. These scripts therefore currently assume that the input
`full_table.tsv` files were generated with the corresponding BUSCO `odb10`
database used to build that reference table.

If you run BUSCO with a different database version, BUSCO identifiers will not
match the reference table correctly. 

## Plot appearance

The plotting step in this repository was reimplemented in Python with
`matplotlib`, as an alternative to the original R plotting scripts from
`lep_busco_painter`.

This was done partly for easier integration into genome notes workflows, but
also to improve figure appearance and control. In particular, the Python
plotting code makes it easier to:

- control font selection and fallback behaviour more cleanly
- place Merian labels for each chromosome to improve interpretation of the plot
- tune figure sizing and spacing for genome note outputs
- split large karyotypes across columns while keeping a consistent Mb scale across panels
- write PNG and SVG outputs directly from the same code path
- customise palette choices while keeping the style reproducible.

The underlying BUSCO-to-Merian interpretation is the same, but the plotting
layer is easier to maintain and better suited to publication-style figure
preparation.

## Files

- `buscopainter.py`: map BUSCO hits to Merian elements and optionally fetch chromosome lengths
- `plot_buscopainter.py`: generate PNG and SVG Merian plots
- `busco_to_merian.sh`: batch wrapper for ToLID lists
- `Merian_elements_full_table.tsv`: reference BUSCO-to-Merian table
- `LICENSE`: MIT license retained for the adapted codebase

## Requirements

- Python 3.10+
- `pandas`
- `matplotlib`
- `requests`
- `biopython`

Install with:

```bash
python3 -m pip install -r requirements.txt
```

## Single assembly workflow

Prepare a plot from a BUSCO result directory:

```bash
mkdir -p output/ilHelArmi9

python3 buscopainter.py \
  --reference_table Merian_elements_full_table.tsv \
  --query_table data/ilHelArmi9/full_table.tsv \
  --prefix output/ilHelArmi9/ \
  --accession GCA_963930815.1

python3 plot_buscopainter.py \
  --file output/ilHelArmi9/all_location.tsv \
  --lengths output/ilHelArmi9/chrom_lengths.tsv \
  --prefix output/ilHelArmi9/ilHelArmi9 \
  --minimum 1 \
  --palette merianbow4 \
  --label-threshold 5 \
  --panel-size 40
```

Outputs:

- `all_location.tsv`
- `chrom_lengths.tsv`
- `*.png`
- `*.svg`

If an assembly has many chromosomes, `plot_buscopainter.py` can split the figure
into columns automatically. `--panel-size` sets the maximum number of
chromosomes per column.

## Plotting modes

These scripts are intentionally used in two different ways.

If you run `buscopainter.py` with `--accession`, chromosome lengths are fetched
from the NCBI Datasets `sequence_reports` endpoint. In that mode, the plot is
chromosome-focused: plotting units are assembled chromosomes, and lengths from
unlocalized scaffolds are added to their parent chromosome. This is the mode
intended for public assemblies and genome note figures. In this mode, the
plotted chromosome names come from the GenBank accessions
returned by NCBI Datasets.

If you do not provide `--accession`, curation assemblies can be plotted against
real scaffold lengths by passing a `.fai` file to `plot_buscopainter.py` with
`--lengths`. In this mode, plotted names come directly from the
sequence identifiers in the BUSCO table and `.fai` index, and scaffold
lengths come from the index.

The plotting step accepts either the NCBI-derived `chrom_lengths.tsv` written
by `buscopainter.py` or a standard `.fai` index. Estimating lengths directly
from BUSCO positions is still supported as a last-resort fallback when no
lengths file is available, but it reflects BUSCO extent rather than true
chromosome or scaffold length.

## Batch workflow

The wrapper expects to be run from a directory containing:

- `tolids`: one ToLID per line
- `tolids_accessions.tsv` or `tolid_accessions.tsv`: two-column tab-separated file with `ToLID<TAB>assembly_accession`

Run:

```bash
bash busco_to_merian.sh
```

The wrapper is environment-variable driven. By default it uses:

- `DATA_ROOT=.` 
- `BUSCO_DIR=${DATA_ROOT}/busco`
- `OUTPUT_DIR=${DATA_ROOT}/merians`

You can override any of these:

```bash
DATA_ROOT=/path/to/project_data \
TOLID_FILE=tolids \
ACCESSION_FILE=tolids_accessions.tsv \
MERIAN_REF=Merian_elements_full_table.tsv \
BUSCO_DIR=/path/to/busco \
OUTPUT_DIR=/path/to/output \
bash busco_to_merian.sh
```

If you prefer, keep the variables in a shell file and source them before
running the batch script:

```bash
source .env
bash busco_to_merian.sh
```

## Example data

Single-panel example for `ilHelArmi9` (`GCA_963930815.1`):

- BUSCO table: [examples/ilHelArmi9_full_table.tsv](examples/ilHelArmi9_full_table.tsv)
- plot output: `examples/ilHelArmi9.png`

This example is plotted using the calculated lengths from the NCBI sequence report, and labelled using chromosome accession numbers.

![Example Merian plot for ilHelArmi9](examples/ilHelArmi9.png)

Multi-panel example for `ilApoPilo2`, a curation in progress:

- BUSCO table: [examples/ilApoPilo2_full_table.tsv](examples/ilApoPilo2_full_table.tsv)
- scaffold index: [examples/ilApoPilo2.fa.fai](examples/ilApoPilo2.fa.fai)
- plot output: `examples/ilApoPilo2.png`

This example uses a `.fai` index so all BUSCO-bearing scaffolds are plotted
against their true lengths. Assemblies with many sequences are split across
columns automatically when the number of plotted chromosomes/scaffolds exceeds
`--panel-size` (default `40`).

![Example Merian plot for ilApoPilo2](examples/ilApoPilo2.png)


## Notes

- `buscopainter.py` keeps both `Complete` and `Duplicated` BUSCO records.
- BUSCO inputs should currently be generated with the `odb10` database version
  compatible with `Merian_elements_full_table.tsv`.
- Chromosome lengths come from NCBI Datasets `sequence_reports`, using the main
  assembled-molecule accession and summing any unlocalized scaffolds assigned to
  that chromosome.
- In NCBI-driven plots, chromosome labels are the assembled-molecule GenBank
  accessions returned by Datasets rather than original assembly header names.
- For curation assemblies, `plot_buscopainter.py` can use a `.fai` index to
  plot assembly sequences against their true scaffold lengths.
- Without a lengths file, the fallback plot includes all BUSCO-bearing
  `query_chr` values from the BUSCO table and estimates their lengths from BUSCO
  extent only.
- The plotting script labels each scaffold/chromosome with Merian elements that
  meet the `--label-threshold`.
- Large chromosome sets can be split across columns with `--panel-size`; panel
  widths shrink or expand to fit their local x-axis range while preserving a
  common Mb scale across the whole figure.
- If Open Sans is available locally it will be used automatically; otherwise
  matplotlib's default sans-serif font is used.


## Attribution

This code adapts ideas, scripts and reference data from
[`charlottewright/lep_busco_painter`](https://github.com/charlottewright/lep_busco_painter).

The `MerianBow4` palette used here is credited to [Arnaud Martin](<https://biology.columbian.gwu.edu/arnaud-martin>).

## Citation

This repository uses the Merian element framework described in:

Wright CJ, et al. 2024. Comparative genomics reveals the dynamics of chromosome evolution in Lepidoptera. *Nature Ecology & Evolution*. doi: [10.1038/s41559-024-02329-4](https://doi.org/10.1038/s41559-024-02329-4), PMID: 38383850; PMCID: PMC11009112
