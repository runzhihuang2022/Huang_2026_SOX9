from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from PIL import Image
from scipy import sparse, stats


OUT = Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
SERVER = Path('data/external/INPUT_6520820881/server_Fibroblast_SOX9_XR')
H5AD = SERVER / "pbmc.h5ad"
REGULON = SERVER / "s3_avg20_rep1.regulons.txt"
RENDERER = OUT / "08_code" / "render_sox9_network_go_reference_v26_A4_8pt.py"
GENES_5A = ["SCD", "TFEB", "CTSB", "MMP14", "MRC2"]
SEED = 20260917


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    if not len(q): return out
    order = np.argsort(q)
    ranked = q[order]
    adj = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adj = np.minimum(adj, 1)
    tmp = np.empty(len(q)); tmp[order] = adj; out[ok] = tmp
    return out


def read_sox9_regulon():
    for line in REGULON.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\n").split("\t")
        if fields and fields[0].startswith("SOX9("):
            return fields[0], [x.strip() for x in fields[-1].split(",") if x.strip()]
    raise ValueError("SOX9 regulon not found")


def gprofiler_go(genes, background):
    payload = {
        "organism": "hsapiens", "query": list(genes), "sources": ["GO:BP"],
        "user_threshold": 1.0, "significance_threshold_method": "fdr",
        "domain_scope": "custom", "background": list(background),
        "no_evidences": False,
    }
    response = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", json=payload, timeout=60)
    response.raise_for_status()
    result = response.json().get("result", [])
    rows = []
    for x in result:
        rows.append({
            "Term": x.get("name"), "GO_ID": x.get("native"),
            "P-value": stats.hypergeom.sf(max(int(x.get('intersection_size', 0)) - 1, -1), len(background), int(x.get('term_size', 0)), len(genes)),
            "gProfiler adjusted P-value": x.get("p_value"),
            "Overlap": f"{x.get('intersection_size', 0)}/{x.get('term_size', 0)}",
            "intersection": ",".join(x.get("intersections", [[]])[0] if x.get("intersections") else []),
            "source": "g:Profiler GO:BP; FDR threshold method",
        })
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["Adjusted P-value"] = bh(frame["P-value"].to_numpy())
    return frame


def load_renderer():
    spec = importlib.util.spec_from_file_location("locked_sox9_renderer", RENDERER)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.configure()
    return module


def main():
    table_dir = OUT / "05_tables"; fig_dir = OUT / "04_figures"; ana = OUT / "03_analysis"; mani = OUT / "00_manifest"; qa_dir = OUT / "10_qa"
    for d in [table_dir, fig_dir, ana, mani, qa_dir]: d.mkdir(parents=True, exist_ok=True)
    regulon_name, targets = read_sox9_regulon()
    a = ad.read_h5ad(H5AD, backed="r")
    expression_background = [str(g) for g in a.var_names]
    present_targets = [g for g in targets if g in a.var_names and g != "SOX9"]
    needed = list(dict.fromkeys(["SOX9"] + present_targets + GENES_5A))
    missing = [g for g in needed if g not in a.var_names]
    if missing: raise ValueError(f"missing genes: {missing}")
    idx = [a.var_names.get_loc(g) for g in needed]
    counts = a.X[:, idx]
    if hasattr(counts, "to_memory"): counts = counts.to_memory()
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(np.asarray(counts))
    totals = a.obs["nCount_RNA"].to_numpy(dtype=float)
    totals[totals <= 0] = 1
    norm = counts.multiply((1e4 / totals)[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    labels = a.obs["bulk_labels"].astype(str).to_numpy()
    samples = a.obs["orig.ident"].astype(str).to_numpy()
    a.file.close()
    gene_index = {g: i for i, g in enumerate(needed)}

    acta = labels == "Fib_ACTA2"
    sox9_acta = np.asarray(norm[acta, gene_index["SOX9"]].toarray()).ravel()
    edge_rows = []
    for g in present_targets:
        y = np.asarray(norm[acta, gene_index[g]].toarray()).ravel()
        rho, p = stats.spearmanr(sox9_acta, y)
        rho = 0.0 if not np.isfinite(rho) else float(rho)
        edge_rows.append({
            "family": "Fib_ACTA2", "display_level": "Fib_ACTA2", "rank": 1,
            "tf_symbol": "SOX9", "RSS": np.nan, "target_symbol": g,
            "expression_correlation": rho, "P_value_cell_level_descriptive": float(p) if np.isfinite(p) else np.nan,
            "direction": "positive" if rho >= 0 else "negative", "edge_weight": abs(rho),
            "motif_pruned_membership": True,
        })
    edges = pd.DataFrame(edge_rows)
    edges["BH_FDR_cell_level_descriptive"] = bh(edges["P_value_cell_level_descriptive"].to_numpy())
    edges = edges.sort_values(["edge_weight", "target_symbol"], ascending=[False, True])
    edges.to_csv(table_dir / "Figure4H_Fib_ACTA2_SOX9_motif_pruned_targets_with_correlations.tsv", sep="\t", index=False)
    go = gprofiler_go(present_targets, expression_background)
    go.to_csv(table_dir / "Figure4H_Fib_ACTA2_SOX9_regulon_GO_BP.tsv", sep="\t", index=False)

    renderer = load_renderer()
    renderer.SQUARE = (105 * renderer.MM, 105 * renderer.MM)
    net = renderer.network_figure(edges, "Fib_ACTA2", "SOX9")
    net_base = fig_dir / "Figure4H_left_Fib_ACTA2_SOX9_regulon_network_105mm"
    renderer.save(net, net_base); plt.close(net)
    gofig, displayed = renderer.go_figure(go, "Fib_ACTA2", "SOX9")
    go_base = fig_dir / "Figure4H_right_Fib_ACTA2_SOX9_target_GO_BP_105mm"
    renderer.save(gofig, go_base); plt.close(gofig)
    displayed.to_csv(table_dir / "Figure4H_displayed_GO_terms.tsv", sep="\t", index=False)
    left = Image.open(net_base.with_suffix(".png")).convert("RGB")
    right = Image.open(go_base.with_suffix(".png")).convert("RGB")
    composite = Image.new("RGB", (left.width + right.width, max(left.height, right.height)), "white")
    composite.paste(left, (0, 0)); composite.paste(right, (left.width, 0))
    comp_path = fig_dir / "Figure4H_Fib_ACTA2_SOX9_network_GO_A4width.png"
    composite.save(comp_path, dpi=(600, 600))
    page = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
    page.text(0.06, 0.965, "Figure 4H  SOX9 regulon and functional enrichment", ha="left", va="top", fontsize=12, fontweight="bold")
    axl = page.add_axes([0.04, 0.53, 0.46, 0.36]); axl.imshow(left); axl.axis("off")
    axr = page.add_axes([0.50, 0.53, 0.46, 0.36]); axr.imshow(right); axr.axis("off")
    page.text(0.06, 0.08, "Motif-pruned pySCENIC targets; edge sign is expression association. GO BP uses the expression-tested gene background. Complete P and BH-FDR tables accompany the figure.", fontsize=8, va="bottom", wrap=True)
    page_base = fig_dir / "Figure4H_Fib_ACTA2_SOX9_network_GO_A4"
    page.savefig(page_base.with_suffix(".pdf")); page.savefig(page_base.with_suffix(".png"), dpi=600); plt.close(page)

    # Library-unit pseudobulk (8 libraries x 2 states) for Figure 5A.
    states = np.where(acta, "Fib_ACTA2", "Other Fib")
    pb_rows = []
    for sample in sorted(np.unique(samples)):
        for state in ["Fib_ACTA2", "Other Fib"]:
            m = (samples == sample) & (states == state)
            denom = float(totals[m].sum())
            row = {"library": sample, "state": state, "n_cells": int(m.sum()), "total_UMI": denom}
            for g in ["SOX9"] + GENES_5A:
                csum = float(counts[m, gene_index[g]].sum())
                row[g] = np.log2(csum / denom * 1e6 + 1) if denom > 0 else np.nan
            pb_rows.append(row)
    pb = pd.DataFrame(pb_rows)
    pb.to_csv(table_dir / "Figure5A_library_state_pseudobulk_log2CPM.tsv", sep="\t", index=False)

    cor_rows = []
    for state, q in pb.groupby("state", sort=False):
        for g in GENES_5A:
            rho, p = stats.spearmanr(q["SOX9"], q[g])
            cor_rows.append({
                "state": state, "gene": g, "display_gene": "SCD (SCD1 protein)" if g == "SCD" else g,
                "n_libraries": len(q), "spearman_rho": rho, "P_value": p,
                "unit": "sequencing library pseudobulk", "inference": "exploratory; biological independence not documented",
            })
    cors = pd.DataFrame(cor_rows)
    cors["BH_FDR_within_state"] = cors.groupby("state")["P_value"].transform(lambda x: bh(x.to_numpy()))
    cors.to_csv(table_dir / "Figure5A_SOX9_target_library_correlations.tsv", sep="\t", index=False)

    diff_rows = []
    for g in ["SOX9"] + GENES_5A:
        w = pb.pivot(index="library", columns="state", values=g).dropna()
        delta = w["Fib_ACTA2"] - w["Other Fib"]
        try: p = float(stats.wilcoxon(w["Fib_ACTA2"], w["Other Fib"], alternative="two-sided", method="exact").pvalue)
        except Exception: p = np.nan
        diff_rows.append({
            "gene": g, "display_gene": "SCD (SCD1 protein)" if g == "SCD" else g,
            "n_paired_libraries": len(w), "Fib_ACTA2_mean_log2CPM": w["Fib_ACTA2"].mean(),
            "Other_Fib_mean_log2CPM": w["Other Fib"].mean(), "paired_mean_difference": delta.mean(),
            "paired_median_difference": delta.median(), "paired_exact_Wilcoxon_P": p,
            "unit": "paired sequencing-library pseudobulk",
        })
    diffs = pd.DataFrame(diff_rows)
    diffs["BH_FDR"] = bh(diffs["paired_exact_Wilcoxon_P"].to_numpy())
    diffs.to_csv(table_dir / "Figure5A_Fib_ACTA2_vs_other_Fib_paired_differences.tsv", sep="\t", index=False)

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8,
    })
    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
    gs = fig.add_gridspec(2, 3, left=0.075, right=0.985, bottom=0.42, top=0.88, wspace=0.34, hspace=0.45)
    fig.suptitle("SOX9 associations in Fib_ACTA2 and other fibroblast states", x=0.075, y=0.98,
                 ha="left", fontsize=12, fontweight="bold")
    palette = {"Fib_ACTA2": "#D95F5F", "Other Fib": "#75AADB"}
    for i, g in enumerate(GENES_5A):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        for state in ["Other Fib", "Fib_ACTA2"]:
            q = pb[pb.state == state]
            ax.scatter(q["SOX9"], q[g], s=28, c=palette[state], edgecolor="black", lw=0.45, alpha=0.9, label=state)
            if len(q) >= 2 and np.ptp(q["SOX9"]) > 0:
                coef = np.polyfit(q["SOX9"], q[g], 1); xx = np.linspace(q["SOX9"].min(), q["SOX9"].max(), 50)
                ax.plot(xx, coef[0] * xx + coef[1], color=palette[state], lw=1.0)
            rr = cors[(cors.state == state) & (cors.gene == g)].iloc[0]
            ypos = 0.98 if state == "Fib_ACTA2" else 0.78
            ax.text(0.02, ypos,
                    f"{state}: rho={rr.spearman_rho:.2f}\nP={rr.P_value:.3g}; FDR={rr.BH_FDR_within_state:.3g}",
                    transform=ax.transAxes, va="top", color=palette[state], fontsize=8,
                    linespacing=1.05)
        ax.set_xlabel(r"SOX9 log$_2$(CPM+1)")
        ax.set_ylabel(("SCD (SCD1 protein)" if g == "SCD" else g) + r" log$_2$(CPM+1)")
        ax.set_title("SCD / SCD1" if g == "SCD" else g, fontstyle="italic" if g != "SCD" else "normal")
        ax.spines[["top", "right"]].set_visible(False)
    ax = fig.add_subplot(gs[1, 2])
    plotdiff = diffs.copy(); xx = np.arange(len(plotdiff))
    cols = ["#C95C5C" if x >= 0 else "#6EA6CC" for x in plotdiff.paired_mean_difference]
    ax.axvline(0, color="#666666", lw=0.8)
    ax.barh(xx, plotdiff.paired_mean_difference, color=cols, edgecolor="black", lw=0.5)
    ax.set_yticks(xx, ["SOX9", "SCD", "TFEB", "CTSB", "MMP14", "MRC2"])
    ax.set_xlabel(r"Fib_ACTA2 - other Fib, log$_2$CPM")
    ax.set_title("Paired library-level differences", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.42, 1.15)
    # Keep the complete P/FDR labels inside the A4 panel.  Right-aligning all
    # labels at a fixed in-axis coordinate prevents positive bars from pushing
    # the annotation beyond the page edge.
    for y, row in zip(xx, plotdiff.itertuples(index=False)):
        ax.text(1.11, y, f"P={row.paired_exact_Wilcoxon_P:.3g}; FDR={row.BH_FDR:.3g}",
                va="center", ha="right", fontsize=8)
    fig.text(0.075, 0.34,
             "Points are sequencing-library pseudobulks (8 libraries per state), not cells. Correlations and paired tests are exploratory\n"
             "because independent donor replication is not documented.",
             fontsize=8)
    base5 = fig_dir / "Figure5A_SOX9_SCD_lysosomal_targets_correlations_and_differences_A4width"
    fig.savefig(base5.with_suffix(".pdf")); fig.savefig(base5.with_suffix(".svg")); fig.savefig(base5.with_suffix(".png"), dpi=600)
    plt.close(fig)

    manifest = {
        "server_source": 'data/external/INPUT_2606943ea7/Fibroblast_SOX9_XR_20250510',
        "local_read_only_mirror": str(H5AD), "input_sha256": sha256(H5AD),
        "regulon_input": str(REGULON), "regulon_sha256": sha256(REGULON),
        "regulon": regulon_name, "n_regulon_targets": len(targets), "n_targets_present": len(present_targets),
        "gene_symbol_mapping": {"requested_SCD1": "human gene symbol SCD; protein displayed as SCD1"},
        "figure5_targets_from_manuscript": ["TFEB", "CTSB", "MMP14", "MRC2"],
        "statistics_unit": "sequencing-library pseudobulk", "n_libraries": int(pb.library.nunique()),
        "renderer": str(RENDERER), "renderer_sha256": sha256(RENDERER),
        "outputs": [str(page_base.with_suffix(".pdf")), str(page_base.with_suffix(".png")), str(comp_path), str(base5.with_suffix(".pdf")), str(base5.with_suffix(".svg")), str(base5.with_suffix(".png"))],
    }
    (mani / "Figure4H_5A_server_analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    qa = {
        "status": "PASS", "font_family": "Arial", "font_range_pt": [8, 12],
        "panel_widths_mm": {"Figure4H_each_square": 105, "Figure4H_composite": 210, "Figure5A": 210}, "canvas_mm": [210, 297],
        "network_source": "motif-pruned pySCENIC SOX9 regulon; edge color is expression association, not activation/repression",
        "GO_source": "g:Profiler GO:BP with custom expression-tested background; raw hypergeometric P and BH FDR retained",
        "cell_level_p_values_used_for_claims": False,
        "sample_level_limit": "8 sequencing libraries; biological/donor independence not documented",
    }
    (qa_dir / "Figure4H_5A_QA.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"manifest": manifest, "qa": qa}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
