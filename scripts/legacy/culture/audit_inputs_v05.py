from pathlib import Path
import hashlib
import json
import anndata as ad
import pandas as pd

OUT = Path('data/external/INPUT_069ac27a1a/01_input_audit')
INPUTS = {
    "human": Path('data/external/INPUT_decd995edd/pbmc_final.h5ad'),
    "rat": Path('data/external/INPUT_a5f4df094a/pbmc_final.raw_counts.h5ad'),
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()

rows = []
details = {}
for species, path in INPUTS.items():
    a = ad.read_h5ad(path, backed="r")
    wanted = [x for x in ["bulk_labels", "bulk_labels.y.y", "sample_id", "sample", "Grade"] if x in a.obs]
    fields = {}
    for key in wanted:
        vc = a.obs[key].astype(str).value_counts(dropna=False)
        fields[key] = {str(k): int(v) for k, v in vc.items()}
    details[species] = {
        "path": str(path), "sha256": sha256(path), "shape": list(a.shape),
        "layers": list(a.layers.keys()), "obsm": list(a.obsm.keys()),
        "obs_columns": list(map(str, a.obs.columns)), "field_counts": fields,
    }
    for key, counts in fields.items():
        for value, n in counts.items():
            rows.append({"species": species, "field": key, "value": value, "n_cells": n})
    a.file.close()

(OUT / "input_metadata_audit.json").write_text(json.dumps(details, ensure_ascii=False, indent=2), encoding="utf-8")
pd.DataFrame(rows).to_csv(OUT / "input_metadata_counts.tsv", sep="\t", index=False)
print(json.dumps(details, ensure_ascii=False, indent=2))
