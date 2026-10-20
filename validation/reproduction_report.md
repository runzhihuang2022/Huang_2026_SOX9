# Reproduction report — 16 October 2026

This report distinguishes numerical replay, complete Fib reclustering, downstream cultured-Fib calculation, and missing upstream/wet-lab steps.

- aggregate: executed_verified; 35 numerical checks; 7.11 s.
- full_culture: executed_verified; 14 numerical checks; 25.07 s.
- full_figure2: executed_verified; 7 numerical checks; 60.58 s.
- source_spatial: executed_verified; 32 numerical checks; 34.24 s.
- Independent-directory aggregate replay: executed_verified; 35 checks; 9 output CSVs byte-identical.
- The clean replay contained only repository files and resolved no external input. The operating system did not sandbox/mount-hide other drives; independence is verified by code review, path scan, and the self-contained executed entry point, not by a filesystem isolation claim.
- Software installation in a brand-new environment and actual GitHub Actions execution were not tested. Tests used the recorded scientific Python environment.
- Source-H5AD spatial gates and human regions were rerun on ten inputs (32 checks). Full raw-data extraction, 465450-cell atlas refit, raw spatial-object map rerendering, native Seurat pipeline, pySCENIC inference, and wet-lab statistics are NOT newly verified.
- Cached GO/network rendering is downstream replay only. Empty KEGG, unresolved term-gene intersections and Venn membership remain open.
- Original manuscript panel letters differ from some later requested computational slots. Use analysis IDs in the coverage ledger until the final legend is reconciled.
- Direct software pins are not complete transitive lockfiles. No third-party font binaries are redistributed.

Full-data outputs containing per-cell expression/identifiers are retained only in the separate author-side private audit directory.

Read `provenance/missing_data.csv` before asserting complete manuscript reproduction. Publication readiness remains partial pending source data, rights/license and author metadata confirmation.
