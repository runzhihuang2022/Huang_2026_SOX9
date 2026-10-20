# Full analysis versus processed-data replay

The default `all` task is intentionally self-contained. The following opt-in
tasks require authorized data obtained separately. Never commit those objects.

## Figure 2A–E: complete lineage clustering

```bash
python -m pip install -r environment/requirements-full.txt
python scripts/reproduce.py full-figure2 --input data/external/selected19_full_QC_Fib_normalized_RNA.h5ad --sample-key public_subject_id --out outputs/full_fib_01
```

Input schema: 20,270 cells × 35,379 genes; source-normalized nonnegative RNA in
`layers['source_normalized_RNA']`; original `bulk_labels == 'Fibroblast'`;
`condition_pair` = normal/fibrosis; pseudonymous `public_subject_id` preserving
the approved patient aggregation. Cell names must be unique and consistently
ordered. If working in a controlled local environment, `--sample-key` may name
the source subject column. Do not publish that column or its key.

The executable in `scripts/renderers/run_figure2_ae_universal.py` is the pristine
project-specific skill executable, copied without analytical or plotting edits.
Its SHA-256 is recorded. It selects 3,000 Seurat HVGs, computes PCA, neighbors,
Leiden and UMAP, freezes geometry, emits A–E panels and marker/sample tables.
Only one group contrast is registered; the unregistered second contrast must
remain NOT_ESTIMABLE. KRT14 is excluded from marker display only, not the RNA.
The SOX9/COL1A1 >0 flag is auxiliary, never the Fib inclusion gate.

Input preparation/QC code is archived in `scripts/legacy/atlas/`. Its absolute
paths have been replaced with external-data aliases. These files preserve the
historical implementation; they are **not a complete standalone raw-data
workflow without the protected metadata and staging map**. The preparation
adapter enforces 200–8,000 features, 500–60,000 counts, ≤20% mitochondrial and
≤1% hemoglobin values, together with the source analysis-QC flag. It retains
208,532 of 208,616 skin cells and then all 20,270 original Fib cells.

## Cross-organ embedding refit

```bash
python scripts/reproduce.py full-atlas --input data/external/atlas_normalized_465450x1679.h5ad --out outputs/full_atlas_01
```

This starts from the approved normalized 1,679-gene analysis matrix, not raw
FASTQ or raw counts. The command clears old embeddings/graphs and repeats the
original centered ARPACK PCA (50), neighbors (20), UMAP (seed 20261009, min_dist
0.35, 200 epochs), and Leiden (0.5). It has not been executed again during this
release unless the validation report explicitly states otherwise.

The archived `build_full_atlas.py` records the upstream source-specific full
count normalization, shared unique uppercase-symbol rule, and feature union.
Original metadata and source matrices are excluded. Uppercase symbol matching
is not validated one-to-one orthology; the embedding is descriptive.

## Cultured Fib matrix-level reconstruction

```bash
python scripts/reproduce.py full-culture --input data/external/culture_Fib_counts.h5ad --out outputs/full_culture_01
```

Input: 59,854 annotated cells × 24,884 genes; nonnegative counts in `X`;
`nCount_RNA` is the full-gene UMI denominator; `bulk_labels` includes Fib_ACTA2;
`orig.ident` records the eight treatment libraries. Source names beginning
`pbmc` are historical object names, not PBMC sample identities.

This runs the archived project's expression normalization, 37-target edge
correlations, pseudobulk calculations, and Figure 4/5 rendering. It freezes the
original cached GO export rather than requesting a changing web database.
GO membership and original pySCENIC upstream execution remain incomplete
evidence, even if the downstream calculations reproduce.

The source RDS has 59,855 entries including one unannotated entry; do not silently
equate that total with the 59,854 annotated H5AD cells. Original native Seurat
analysis code/lockfile was not located in the inspected source directory.
Observed R package versions are not retrospective proof of the original run.

## Spatial source-level reconstruction

The source reading/gating/rendering implementations are archived in
`scripts/legacy/spatial/` and `scripts/renderers/spatial_primitives.py`.
`data/source_locations.csv` lists expected hashes and access status; author-held
staging information is intentionally outside this repository. Source X is used
as stored, with no second normalization. Each gene must exist exactly once;
the parent Fib annotation is preserved; double positive means both values >0.

Use the exact registered human scar/scarless two-point anchors after the stored
rotation. Reuse the original ROI; never select a new zoom to increase positives.
Rat rendering uses the recorded coordinate key (usually X_rot_spatial) and its
registered ROI, with no additional rotation. Physical scale is the historical
0.5 µm per coordinate unit; it is recorded, not newly calibrated.

The public default validates section-level/region-level counts, not raw H5AD
coordinates. Source-level marker gating and human region assignments can be run:

```bash
python scripts/full_spatial.py --input-dir data/external/spatial --out outputs/full_spatial_01
```

Stage the ten authorized H5AD files using filenames and SHA-256 values in
`config/spatial_sections.json`. That config contains only the registered anchor
geometry, source hashes, coordinate keys and available ROIs, not human barcodes
or RNA. Alternatively, `--sources-json` can read an author-held private mapping
of section IDs to existing source paths; never commit the mapping.

The release preparation executed this source-level command on all ten inputs
and verified 32 checks. Add `--maps` to render source-part full/ROI views with
the approved spatial primitives. That optional map mode was not executed during
release preparation and is not a replacement for the approved manuscript canvas.
Some human source parts have only a registered pooled-display ROI, not an
individual-part ROI; no new zoom is invented for those parts. Full reconstruction
of the original pooled-display canvas still needs its registered geometry parent.
Raw Stereo-seq alignment/cell segmentation/QC is not recreated by this command.
Approved main map PDFs/PNGs are reference outputs, not falsely labeled fresh maps.

## Wet-lab panels and supplementary figures

Raw qPCR values/Ct, uncropped immunoblots and densitometry, growth curves,
migration/contraction measurements, animal IDs/assignments, microscopy counting
rules, histology measurements, and the corresponding Prism/R scripts were not
found in the currently inspected analysis directories. No such values were
digitized from the composite TIFFs. See the missing-data ledger.
