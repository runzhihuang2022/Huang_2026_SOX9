# Figure 1

The global embedding is frozen for default display replay. Welch/BH expression statistics are recomputed. Human main maps are reference outputs; source-level gating has a separate full-data entry point.

## Run

```bash
python Figure1/run_figure.py --out outputs/my_figure1_01
```

Run from repository root. The output directory must be new. Full requirements and data staging are in [FULL_ANALYSIS](../docs/FULL_ANALYSIS.md).

## Processed source data

[Source tables](source_data/) are exact copies of the corresponding canonical aggregate inputs under [data/processed](../data/processed/). They are not new analyses or raw patient data. Shared analysis code is under [scripts](../scripts/).

## Reference panel masters

- [Figure1bc_horizontal_selected19_fullcells.pdf](outputs/Figure1bc_horizontal_selected19_fullcells.pdf)
- [Figure1E_COLI_SOX9_unified_slot.pdf](outputs/Figure1E_COLI_SOX9_unified_slot.pdf)

Panel identifiers follow the current requested computational slots. Final manuscript legend numbering still needs reconciliation with the older TIFFs. Wet-lab panels are not certified by these computational outputs.
