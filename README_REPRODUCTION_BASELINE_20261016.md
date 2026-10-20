> Historical local-preparation snapshot (16 October 2026). For current publication status, use [README.md](README.md).

# SOX9 skin repair and cross-organ fibrosis — reproducibility materials

Local release candidate v01, prepared 16 October 2026. This directory is the
repository root. No repository has been created, pushed, or published.

This package records the **current 9 October analysis**, not the superseded
8 October figure TIFFs. It contains executable analysis/plotting code,
deidentified aggregate inputs, reference computational panels, parameter files,
environment records, checksums, and an explicit coverage/gap ledger.

## Quick start

Use Python 3.10 in a separate environment. From this directory:

```bash
python -m pip install -r environment/requirements-replay.txt
python scripts/validate_release.py
python scripts/reproduce.py all --out outputs/my_replay_01
```

The output directory **must not already exist**. The default command is offline:
it makes no web, SSH, GitHub, or Office calls and does not modify source inputs.
On the preparation workstation the numerical replay takes approximately 10 seconds;
package installation, full clustering, and raw-data extraction are separate tasks.

For the original analysis environment and opt-in matrix-level execution, see
[full-analysis instructions](docs/FULL_ANALYSIS.md). A small numerical replay is
not equivalent to verifying every raw-data-to-manuscript step.

## Frozen scientific scope

| Analysis | Current scope |
|---|---|
| Figure 1b–c | Five organs; 465,450 cells from 96 libraries; no PBMC specimens |
| Skin in Figure 1 and Figure 2 | The same 19 libraries; 208,532 QC-pass skin cells; all 20,270 source-annotated Fib cells |
| Skin expression endpoint | Mean normalized *SOX9* within **all Fib**, including genuine zeros; not SOX9-positive percentage |
| Skin statistical units | 5 normal and 10 recorded scar patient units; equal library weighting within recorded patient |
| Figure 2 clustering | 17 unsupervised C labels; C12 is the highest observed mean *SOX9* cluster |
| C12 detection | 61/185 SOX9-positive (32.973%); ACTA2 180/185; VIM 178/185 |
| Human spatial | 2 mph and pooled-display 6 mph parts in the main panel; 9 mph is supplemental checking |
| Rat spatial | Adult 1/5/19 dpw and E16.5 6 h/3/6 dpw; six sections, not six verified independent animals |
| Cultured Fib | OE1/2, OEC-1/2, SI905-1/2, SINC-1/2; eight treatment libraries |
| Table S1 | 104 actual scRNA library records: 96 atlas + 8 cultured Fib; Figure 2 libraries are not added twice |

The locked 19 skin record IDs are in `config/skin_19_record_whitelist.csv`.
Scalp references and historical excluded libraries are not included in this cohort.
The source parent object's historical name contains `pbmc`; that filename is
**not evidence that PBMC samples were used**. The released sample summary records
only specimens actually used in the current analysis.

## What is reproduced

`reproduce.py all` recomputes the five Welch tests and their joint BH adjustment,
checks the C12 detection percentages, executes the original cluster-composition
and paired-cluster test functions, recomputes section/region counts from aggregate
tables, and recalculates cultured-Fib library-level correlations and paired tests.
It rerenders the original Figure 1 plotting block and Figure 5 plotting block,
the pristine Figure 2 dotplot body, and the existing SOX9 network/GO renderer.

Public UMAP display textures are **frozen analysis-derived rasters**, generated
directly from the approved geometry. They are not digitized experimental images
or original measurements. The default run does not refit UMAP from those rasters.
Per-cell human RNA matrices, cell barcodes/coordinates, and patient linkage keys
are excluded. Full clustering/embedding uses separately authorized inputs.
The Figure 2 aggregate replay uses a frozen UMAP texture, so it is not a
pixel-identical replacement for the approved main-panel master. Its dotplot
formula and numerical inputs are unchanged. Approved masters are under
`outputs/reference/`.

Spatial main maps are supplied as approved reference panels. Their counts are
retested offline. Source-object marker gating and human region assignments were
also rerun on all ten registered human/rat inputs (32 checks). Raw-object spatial
map regeneration requires the external matrices and registered anchors.
Existing experimental panels are not silently
treated as reproduced merely because an old composite TIFF exists.

## Important evidence status

- C12 contains 183 source-annotated Myo_Fib cells and 2 Inflama_Fib cells, but is
  KRT14-high. Removing KRT14 from the dotplot did not remove it from the RNA or
  resolve lineage purity. Do not label C12 a proven pure myofibroblast population.
- The selected skin cohort and post-ranking cluster comparisons are exploratory.
  Cells are not biological replicates. Absence of a cell population is not zero
  gene expression; genuine zero-expression observations are retained.
- Human scar/scarless sides use registered anatomical anchors. Rat wound-centre
  versus flank annotations are not automatically scar-inside/outside masks.
- Spatial comparisons lack independent donor/animal mapping and matched times;
  inferential P and FDR are `NE`. No significant embryo-versus-adult claim is certified.
- Network edges use an existing motif-pruned target list and expression
  association; the original GRNBoost2 → cisTarget → AUCell execution chain is
  not completely archived here. Edge sign is not proven activation/repression.
- The cached regulon GO export has a problematic `intersection` field (some
  entries are evidence codes, not gene symbols); do not treat it as validated
  term-membership evidence. It is preserved and flagged, not silently repaired.
- The source `Fib_ACTA2_KEGG_Enrichment_Up.csv` is empty (3 bytes). KEGG enrichment
  and the manuscript's 39-gene Venn intersection are not certified by this package.
- Original Seurat marker `p_val_adj` values are **Bonferroni**, not BH FDR.
  Source R versions are runtime observations, not an original analysis lockfile.
- The legacy Figure 5A reference PDF contains non-Arial mathematical glyphs.
  It is retained unchanged and is not certified as an all-Arial final figure;
  the per-file font audit records this exception. Numerical reproduction does
  not certify final publication typography.

See [coverage](provenance/figure_manifest.csv),
[missing evidence](provenance/missing_data.csv), and
[executed reproduction report](validation/reproduction_report.md).
Figures 6–7 and wet-lab panels require original experiment-level measurements,
sample assignments, and analysis files before their statistics can be reproduced.
Supplements are inventoried; unlocated analyses are explicitly marked missing.

## Directory map

```text
config/               locked cohort, parameters, full-analysis config
data/processed/       small aggregate numerical source tables
data/metadata/        public sample summary and data dictionary
data/frozen_maps/     approved geometry-derived UMAP display textures
environment/          Python pins and observed R versions
scripts/reproduce.py  documented portable entry point
scripts/renderers/    original configurable renderer / approved routines
scripts/legacy/       path-redacted source archive; not all directly executable
outputs/reference/    approved computational panels (not full experimental TIFFs)
provenance/           coverage, checksums, code manifest, missing-data ledger
validation/           actual test results, privacy/path/size audit
docs/                 full-run staging, scientific boundaries, upload instructions
```

## Data access, privacy, rights, and citation

The public sample summary removes source clinical library names, original patient
IDs, patient linkage keys, sex, anatomical detail, and individual sampling times.
Anonymous statistical-unit labels retain only the grouping needed for the
reported aggregate calculations. No human raw matrices, server addresses,
credentials, private keys, or local access routes are included.

Public source accessions: GSE135893, GSE196929, GSE132146, GSE153480, GSE136103.
Use the original deposit and its terms. Study-generated matrices and full clinical
metadata require the corresponding author's authorization; a public accession
has not been invented. Removing identifiers is not proof of consent for public
release. Authors must confirm aggregate-data and figure reuse rights before upload.

No software/data license, author list, DOI, or public repository URL is assigned
on the authors' behalf. See `LICENSE.proposed.txt` and `CITATION.proposed.md`.
Do not claim this local candidate is already publicly available or that it
reproduces every panel of the manuscript.
