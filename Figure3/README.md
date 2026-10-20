# Figure 3

The runner recomputes six-section/region count validation, not animal-level inference. Existing spatial main maps are approved reference outputs. Source-H5AD gating is available in scripts/full_spatial.py.

## Run

```bash
python Figure3/run_figure.py --out outputs/my_figure3_01
```

Run from repository root. The output directory must be new. Full requirements and data staging are in [FULL_ANALYSIS](../docs/FULL_ANALYSIS.md).

## Processed source data

[Source tables](source_data/) are exact copies of the corresponding canonical aggregate inputs under [data/processed](../data/processed/). They are not new analyses or raw patient data. Shared analysis code is under [scripts](../scripts/).

## Reference panel masters

- [Figure3G_COLI_SOX9_unified_slot.pdf](outputs/Figure3G_COLI_SOX9_unified_slot.pdf)

Panel identifiers follow the current requested computational slots. Final manuscript legend numbering still needs reconciliation with the older TIFFs. Wet-lab panels are not certified by these computational outputs.
