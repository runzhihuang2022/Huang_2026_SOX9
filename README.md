# Huang_2026_SOX9

Analysis code and deidentified processed source data for the SOX9 skin-repair
and cross-organ fibrosis study.

Repository: https://github.com/runzhihuang2022/Huang_2026_SOX9  
Public repository organization prepared: **20 October 2026**  
Scientific baseline: **9 October 2026, selected-19-library full-cell analysis**

## Repository organization

The figure-oriented navigation follows the organizational style of
[Zhang2019NBT](https://github.com/YongxinLiu/Zhang2019NBT) and
[Huang_20251224_iMeta](https://github.com/runzhihuang2022/Huang_20251224_iMeta).
No code, experimental data or licensing terms were copied from those repositories.

| Directory | Contents and scope |
|---|---|
| [Figure1](Figure1) | Cross-organ UMAP reference, SOX9 expression comparisons and human SOX9/COL1A1 spatial panel; processed source tables and runner |
| [Figure2](Figure2) | Same 19 skin libraries; all 20,270 original Fib cells; unsupervised clustering, marker dotplot and sample-level composition analyses |
| [Figure3](Figure3) | Adult/E16.5 rat spatial panel and six-section double-positive count validation |
| [Figure4](Figure4) | Cultured-Fib SOX9 target-network/cached GO material; original native Seurat/KEGG/Venn gaps explicitly recorded |
| [Figure5](Figure5) | Eight-library SOX9/SCD and lysosomal-target association and paired-difference analysis |
| [Figure6](Figure6), [Figure7](Figure7) | Coverage records for experimental panels; raw quantitative sources are not included or claimed reproduced |
| [data](data) | Deidentified aggregate numerical inputs, metadata dictionary, frozen geometry-derived display maps and external-input contracts |
| [script](script), [scripts](scripts) | Convenience launcher and canonical reproducibility/plotting scripts |
| [table](table) | Privacy-minimized Table S1 for the 104 scRNA libraries actually used |
| [config](config), [environment](environment) | Locked selection, analysis parameters, direct software pins and observed runtime versions |
| [outputs/reference](outputs/reference) | Existing approved computational-panel masters, not source experimental measurements |
| [provenance](provenance), [validation](validation) | Source hashes, figure coverage, executed tests and unresolved evidence |

Each implemented figure folder contains a figure-specific runner, small source-data
tables and reference outputs. Shared canonical data/scripts stay in their original
relative-path locations; they are not reimplemented in another language.

## Quick start

```bash
git clone https://github.com/runzhihuang2022/Huang_2026_SOX9.git
cd Huang_2026_SOX9
python -m pip install -r environment/requirements-replay.txt
python scripts/validate_release.py
python scripts/reproduce.py all --out outputs/my_replay_01
```

Use Python 3.10 and licensed Arial for the figure-rendering command. The output
directory must be new. On systems without Arial, numerical validation can run:

```bash
python scripts/reproduce.py all --out outputs/my_tables_01 --tables-only
```

An individual figure can be replayed, for example:

```bash
python Figure1/run_figure.py --out outputs/my_figure1_01
```

Full lineage clustering, count-matrix cultured-Fib calculations and source-H5AD
spatial gating require separately authorized inputs. See
[FULL_ANALYSIS.md](docs/FULL_ANALYSIS.md).

## Verified analysis scope

- Figure 1 contains 465,450 cells from 96 libraries. Its skin cohort contains
  exactly 19 selected libraries and 208,532 QC-pass cells.
- Figure 2 uses the same 19 skin libraries and all 20,270 original source-annotated
  Fib cells. C12 has the highest observed mean SOX9 expression; 61/185 cells
  detect SOX9, 180/185 detect ACTA2, and 178/185 detect VIM.
- Current Table S1 contains 104 actual scRNA libraries: 96 atlas libraries plus
  eight cultured-Fib treatment libraries. PBMC specimens are not included.
- The spatial gate is source Fib annotation AND SOX9 >0 AND COL1A1 >0
  (rat: Sox9 >0 AND Col1a1 >0). Ten human/rat source inputs were checked;
  the main rat panel uses six sections.
- The archived validation records 35 aggregate replay checks, 7 full Fib
  reclustering checks, 14 cultured-Fib matrix-level checks and 32 source-spatial
  checks. Nine replayed output CSVs matched byte-for-byte in an independent
  directory on the recorded environment.

These tests are not certification that all manuscript experiments or all
upstream processing can be reconstructed from the public bundle. See the
[executed report](validation/reproduction_report.md),
[coverage ledger](provenance/figure_manifest.csv) and
[missing-data ledger](provenance/missing_data.csv).

## Data availability and privacy

This repository provides **analysis code and deidentified processed source data**
for the computational analyses documented here. Original human expression
matrices, individual cell barcodes/coordinate matrices, original patient/library
linkage keys, server locations, credentials and private audit records are excluded.
They must not be inferred to be publicly deposited merely because this repository
is public. Study-generated raw-data access/deposition must be stated separately
by the corresponding author.

Public-source accessions used by the current atlas are GSE135893, GSE196929,
GSE132146, GSE153480 and GSE136103. Source deposits retain their original terms.

## Interpretation and completeness

The skin endpoint in Figure 1 is mean normalized SOX9 in **all Fib**, including
true zero expression, not SOX9-positive proportion. Skin tests use recorded
patient-unit means; the selected cohort and post-ranking comparisons are
exploratory. Missing target-cell populations are not coded as expression zero.

C12 is KRT14-high despite its predominantly original Myo_Fib annotation; it must
not be described as a proven pure myofibroblast population. KRT14 was removed
from the dotplot display only, not from the source expression data.

Spatial sections are not verified independent donors/animals and their sampled
times differ. Their P/FDR are NE; this repository does not certify a significant
embryo-versus-adult group effect. The main UMAP/spatial reference panels include
frozen displays; default replay does not refit missing raw matrices from images.

Original Seurat marker p_val_adj is Bonferroni-adjusted, not BH FDR. The source
KEGG-up export is empty; the Venn gene intersection and the original full native
Seurat/pySCENIC execution chain remain incomplete. Cached GO term-gene membership
and legacy Figure5 mathematical font exceptions are documented, not silently
certified. Missing wet-lab measurements are not digitized from composite TIFFs.

## Citation and rights

The final manuscript citation/DOI and a blanket code/data license have not been
assigned on the authors' behalf. See [CITATION.proposed.md](CITATION.proposed.md)
and [LICENSE.proposed.txt](LICENSE.proposed.txt). Third-party dataset, database,
software and figure rights remain with their respective owners.

## Checksums and historical snapshot

[release_checksums.csv](provenance/release_checksums.csv) describes exact tracked
release-file bytes. `.gitattributes` prevents implicit line-ending conversion.
The [16 October preparation README](README_REPRODUCTION_BASELINE_20261016.md)
is retained as a historical snapshot, not the current publication status.
