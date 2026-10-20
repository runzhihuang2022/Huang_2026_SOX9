# Figure 2

The default runner reuses a frozen UMAP texture and the pristine marker-dotplot body; complete 20,270-cell reclustering is separately runnable with authorized input. C12 lineage purity remains unresolved.

## Run

```bash
python Figure2/run_figure.py --out outputs/my_figure2_01
```

Run from repository root. The output directory must be new. Full requirements and data staging are in [FULL_ANALYSIS](../docs/FULL_ANALYSIS.md).

## Processed source data

[Source tables](source_data/) are exact copies of the corresponding canonical aggregate inputs under [data/processed](../data/processed/). They are not new analyses or raw patient data. Shared analysis code is under [scripts](../scripts/).

## Reference panel masters

- [Figure2AB_human_Fib_slot.pdf](outputs/Figure2AB_human_Fib_slot.pdf)

Panel identifiers follow the current requested computational slots. Final manuscript legend numbering still needs reconciliation with the older TIFFs. Wet-lab panels are not certified by these computational outputs.
