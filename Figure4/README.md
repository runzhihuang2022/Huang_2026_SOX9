# Figure 4

The shared cultured-Fib runner renders the SOX9 target network/cached GO and validates library-level calculations. It does not reconstruct the unavailable original native Seurat, full SCENIC, empty KEGG or Venn workflows.

## Run

```bash
python Figure4/run_figure.py --out outputs/my_figure4_01
```

Run from repository root. The output directory must be new. Full requirements and data staging are in [FULL_ANALYSIS](../docs/FULL_ANALYSIS.md).

## Processed source data

[Source tables](source_data/) are exact copies of the corresponding canonical aggregate inputs under [data/processed](../data/processed/). They are not new analyses or raw patient data. Shared analysis code is under [scripts](../scripts/).

## Reference panel masters

- [Figure4H_left_Fib_ACTA2_SOX9_regulon_network_105mm.pdf](outputs/Figure4H_left_Fib_ACTA2_SOX9_regulon_network_105mm.pdf)
- [Figure4H_right_Fib_ACTA2_SOX9_target_GO_BP_105mm.pdf](outputs/Figure4H_right_Fib_ACTA2_SOX9_target_GO_BP_105mm.pdf)

Panel identifiers follow the current requested computational slots. Final manuscript legend numbering still needs reconciliation with the older TIFFs. Wet-lab panels are not certified by these computational outputs.
