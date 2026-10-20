#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import json
import platform
import shutil
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


REQUIRED_DIRS = [
    "00_manifest", "01_input_audit", "02_shared_geometry", "03_analysis",
    "04_figures", "05_tables", "06_html", "07_reproduction", "08_code",
    "09_logs", "10_qa",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty(len(pv), dtype=float)
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return out


def dense_vector(matrix, column: int) -> np.ndarray:
    value = matrix[:, column]
    if sparse.issparse(value):
        return np.asarray(value.toarray()).ravel()
    return np.asarray(value).ravel()


def validate_config(cfg: dict) -> None:
    required = ["input_h5ad", "output_dir", "dataset_id", "lineage_gate", "expression", "embedding", "markers"]
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing config keys: {missing}")
    markers = cfg["markers"]
    for key in ["mandatory", "allowlist", "target_gene"]:
        if key not in markers:
            raise ValueError(f"markers.{key} is required")
    if cfg.get("panel_mode", "target_cluster") not in {"target_cluster", "composition_contrasts"}:
        raise ValueError("panel_mode must be target_cluster or composition_contrasts")
    if cfg.get("panel_mode") == "composition_contrasts" and not cfg.get("group_key"):
        raise ValueError("composition_contrasts requires group_key")


def prepare_output(path: Path, allow_existing_empty: bool) -> None:
    if path.exists():
        contents = list(path.iterdir())
        if contents or not allow_existing_empty:
            raise FileExistsError(f"Refusing to overwrite output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_DIRS:
        (path / name).mkdir(exist_ok=True)


def lineage_mask(obs: pd.DataFrame, gate: dict) -> np.ndarray:
    column = gate.get("column")
    if not column:
        return np.ones(len(obs), dtype=bool)
    if column not in obs:
        raise KeyError(f"Lineage gate column not found: {column}")
    values = obs[column].astype(str)
    mask = np.zeros(len(obs), dtype=bool)
    include = [str(x) for x in gate.get("include", [])]
    prefixes = [str(x) for x in gate.get("prefix", [])]
    if include:
        mask |= values.isin(include).to_numpy()
    for prefix in prefixes:
        mask |= values.str.startswith(prefix).to_numpy()
    if not include and not prefixes:
        raise ValueError("lineage_gate must contain include or prefix, or set column to null")
    return mask


def expression_matrix(adata: ad.AnnData, cfg: dict):
    normalized = cfg.get("normalized_layer")
    counts_layer = cfg.get("counts_layer")
    if normalized:
        if normalized not in adata.layers:
            raise KeyError(f"Normalized layer not found: {normalized}")
        matrix = adata.layers[normalized].copy()
        source = f"layers/{normalized}"
    else:
        if counts_layer:
            if counts_layer not in adata.layers:
                raise KeyError(f"Counts layer not found: {counts_layer}")
            matrix = adata.layers[counts_layer].copy()
            source = f"layers/{counts_layer}->normalize_total/log1p"
        else:
            matrix = adata.X.copy()
            source = "X->normalize_total/log1p"
        temporary = ad.AnnData(X=matrix, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
        sc.pp.normalize_total(temporary, target_sum=float(cfg.get("target_sum", 10000)))
        sc.pp.log1p(temporary)
        matrix = temporary.X
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    return matrix, source


def marker_summary(matrix, labels, genes, var_names) -> pd.DataFrame:
    rows = []
    positions = pd.Index(var_names).get_indexer(genes)
    for cluster in sorted(pd.unique(labels), key=lambda x: int(str(x)[1:])):
        cell_mask = labels == cluster
        for gene, position in zip(genes, positions):
            values = dense_vector(matrix[cell_mask], int(position))
            rows.append({
                "cluster": cluster,
                "gene": gene,
                "mean_expression": float(np.mean(values)),
                "percent_expressed": float(np.mean(values > 0) * 100),
                "n_cells": int(cell_mask.sum()),
            })
    return pd.DataFrame(rows)


def setup_style(style: dict) -> None:
    minimum = float(style.get("font_min_pt", 8))
    maximum = float(style.get("font_max_pt", 12))
    if minimum < 7 or maximum > 12 or minimum > maximum:
        raise ValueError("Visible font range must remain within 7-12 pt and min <= max")
    mpl.rcParams.update({
        "font.family": style.get("font_family", "Arial"),
        "font.size": minimum,
        "axes.titlesize": min(maximum, 10),
        "axes.labelsize": min(maximum, 9),
        "xtick.labelsize": minimum,
        "ytick.labelsize": minimum,
        "legend.fontsize": minimum,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.8,
        "axes.unicode_minus": False,
    })


def clean_umap(ax) -> None:
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_cluster_umap(ax, coordinates, labels, palette, title, highlight=None) -> None:
    if highlight is None:
        for cluster in sorted(palette, key=lambda x: int(x[1:])):
            mask = labels == cluster
            ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=1.4, color=palette[cluster], linewidths=0,
                       rasterized=True, label=cluster)
    else:
        background = labels != highlight
        ax.scatter(coordinates[background, 0], coordinates[background, 1], s=1.2, color="#E4E4E4", alpha=0.65,
                   linewidths=0, rasterized=True)
        mask = labels == highlight
        ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=2.0, color=palette[highlight], linewidths=0,
                   rasterized=True)
    ax.set_title(title, loc="left", fontweight="bold")
    clean_umap(ax)


def sample_cluster_table(samples, groups, labels, target, positive) -> pd.DataFrame:
    source = pd.DataFrame({"sample": samples, "group": groups, "cluster": labels, "target": target, "positive": positive})
    rows = []
    for (sample, group, cluster), block in source.groupby(["sample", "group", "cluster"], sort=False):
        rows.append({
            "sample": sample, "group": group, "cluster": cluster, "n_cells": len(block),
            "mean_target_expression": float(block.target.mean()),
            "positive_fraction": float(block.positive.mean()),
        })
    return pd.DataFrame(rows)


def highest_cluster_tests(summary: pd.DataFrame, highest: str) -> pd.DataFrame:
    rows = []
    high = summary[summary.cluster.eq(highest)].set_index("sample")
    clusters = sorted(summary.cluster.unique(), key=lambda x: int(x[1:]))
    for cluster in [x for x in clusters if x != highest]:
        other = summary[summary.cluster.eq(cluster)].set_index("sample")
        common = high.index.intersection(other.index)
        delta = high.loc[common, "mean_target_expression"] - other.loc[common, "mean_target_expression"]
        p = np.nan
        if len(delta) >= 3 and np.any(delta != 0):
            p = float(stats.wilcoxon(delta, alternative="greater", method="auto").pvalue)
        rows.append({
            "highest_cluster": highest, "comparison_cluster": cluster, "n_paired_samples": len(common),
            "mean_paired_difference": float(delta.mean()) if len(delta) else np.nan,
            "one_sided_paired_Wilcoxon_P": p,
            "status": "ESTIMABLE" if np.isfinite(p) else "DESCRIPTIVE_OR_NOT_ESTIMABLE",
        })
    frame = pd.DataFrame(rows)
    frame["BH_FDR"] = bh(frame["one_sided_paired_Wilcoxon_P"])
    return frame


def contrast_table(summary: pd.DataFrame, comparison: dict) -> pd.DataFrame:
    numerator = str(comparison["numerator"])
    denominator = str(comparison["denominator"])
    clusters = sorted(summary.cluster.unique(), key=lambda x: int(x[1:]))
    # A missing sample x cluster row means zero cells, not a missing sample.
    # Complete the grid before proportions/tests so every sample contributes to
    # every cluster and group-level compositions sum to exactly 100%.
    registry = summary[["sample", "group"]].drop_duplicates()
    grid = registry.assign(_key=1).merge(pd.DataFrame({"cluster": clusters, "_key": 1}), on="_key").drop(columns="_key")
    counts = grid.merge(summary[["sample", "group", "cluster", "n_cells"]],
                        on=["sample", "group", "cluster"], how="left")
    counts["n_cells"] = counts["n_cells"].fillna(0).astype(int)
    totals = counts.groupby(["sample", "group"], as_index=False).n_cells.sum().rename(columns={"n_cells": "total"})
    fractions = counts.merge(totals, on=["sample", "group"], how="left")
    fractions["fraction"] = fractions.n_cells / fractions.total
    rows = []
    for cluster in clusters:
        a = fractions[(fractions.group == numerator) & (fractions.cluster == cluster)].fraction.to_numpy()
        b = fractions[(fractions.group == denominator) & (fractions.cluster == cluster)].fraction.to_numpy()
        p = np.nan
        if len(a) >= 3 and len(b) >= 3:
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        rows.append({
            "comparison": comparison.get("name", f"{numerator}_vs_{denominator}"), "cluster": cluster,
            "numerator": numerator, "denominator": denominator, "n_numerator_samples": len(a),
            "n_denominator_samples": len(b), "mean_fraction_difference": float(np.mean(a) - np.mean(b)) if len(a) and len(b) else np.nan,
            "Mann_Whitney_P": p, "status": "ESTIMABLE" if np.isfinite(p) else "NOT_ESTIMABLE",
        })
    frame = pd.DataFrame(rows)
    frame["BH_FDR"] = bh(frame.Mann_Whitney_P)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description="Universal Figure 2A-E lineage subtype analysis")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    validate_config(cfg)

    input_path = Path(cfg["input_h5ad"])
    output = Path(cfg["output_dir"])
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    prepare_output(output, bool(cfg.get("allow_existing_empty_output", False)))
    setup_style(cfg.get("style", {}))

    source = ad.read_h5ad(input_path)
    mask = lineage_mask(source.obs, cfg["lineage_gate"])
    if int(mask.sum()) < 20:
        raise ValueError(f"Lineage gate retained only {int(mask.sum())} cells; at least 20 required")
    lineage = source[mask].copy()
    matrix, expression_source = expression_matrix(lineage, cfg["expression"])
    lineage.X = matrix

    markers = cfg["markers"]
    mandatory = list(dict.fromkeys(map(str, markers["mandatory"])))
    allowlist = list(dict.fromkeys(map(str, markers["allowlist"])))
    target_gene = str(markers["target_gene"])
    coexpression = list(dict.fromkeys(map(str, markers.get("coexpression", [target_gene]))))
    required_genes = list(dict.fromkeys(mandatory + [target_gene] + coexpression))
    missing = [gene for gene in required_genes if gene not in lineage.var_names]
    if missing:
        raise ValueError(f"Required genes absent from input: {missing}")
    available_allowlist = [gene for gene in allowlist if gene in lineage.var_names]
    if len(available_allowlist) < 3:
        raise ValueError("Fewer than three allowlist genes are present")

    embedding_cfg = cfg["embedding"]
    seed = int(embedding_cfg.get("seed", 0))
    representation = next((key for key in embedding_cfg.get("use_rep_candidates", []) if key in lineage.obsm), None)
    if representation is None:
        sc.pp.highly_variable_genes(lineage, n_top_genes=min(3000, lineage.n_vars), flavor="seurat")
        sc.pp.pca(lineage, n_comps=min(50, lineage.n_obs - 1, lineage.n_vars - 1), use_highly_variable=True,
                  random_state=seed)
        representation = "X_pca"
    sc.pp.neighbors(lineage, use_rep=representation, n_neighbors=int(embedding_cfg.get("n_neighbors", 20)), random_state=seed)
    sc.tl.leiden(lineage, resolution=float(embedding_cfg.get("resolution", 0.55)), key_added="leiden_universal",
                 random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(lineage, random_state=seed, min_dist=float(embedding_cfg.get("min_dist", 0.35)), spread=1.0)
    raw_labels = lineage.obs.leiden_universal.astype(str).to_numpy()
    raw_order = sorted(pd.unique(raw_labels), key=lambda x: int(x))
    rename = {old: f"C{index + 1}" for index, old in enumerate(raw_order)}
    labels = np.array([rename[value] for value in raw_labels], dtype=object)
    coordinates = np.asarray(lineage.obsm["X_umap"], dtype=float)
    cluster_order = sorted(pd.unique(labels), key=lambda x: int(x[1:]))

    marker = marker_summary(matrix, labels, available_allowlist, lineage.var_names)
    selected_rows = []
    for cluster, block in marker.groupby("cluster", sort=False):
        selected = block[~block.gene.isin(mandatory)].sort_values(
            ["mean_expression", "percent_expressed", "gene"], ascending=[False, False, True]
        ).head(3)
        if len(selected) != 3:
            raise ValueError(f"Cluster {cluster} has fewer than three eligible non-mandatory allowlist markers")
        for rank, row in enumerate(selected.itertuples(index=False), 1):
            selected_rows.append({
                "label": cluster, "gene": row.gene, "rank": rank, "mean_expression": row.mean_expression,
                "percent_expressed": row.percent_expressed, "selection_mode": "all_allowlist_top3",
            })
    ledger = pd.DataFrame(selected_rows)
    displayed_genes = list(dict.fromkeys(mandatory + ledger.gene.tolist()))

    target_index = lineage.var_names.get_loc(target_gene)
    target = dense_vector(matrix, target_index)
    target_summary = marker_summary(matrix, labels, [target_gene], lineage.var_names)
    highest = str(target_summary.sort_values(["mean_expression", "percent_expressed"], ascending=False).iloc[0].cluster)
    co_positions = pd.Index(lineage.var_names).get_indexer(coexpression)
    co_matrix = matrix[:, co_positions]
    threshold = float(markers.get("positive_threshold", 0))
    if sparse.issparse(co_matrix):
        positive = np.asarray((co_matrix > threshold).sum(axis=1)).ravel() == len(co_positions)
    else:
        positive = (np.asarray(co_matrix) > threshold).all(axis=1)

    sample_key = cfg.get("sample_key")
    group_key = cfg.get("group_key")
    samples = lineage.obs[sample_key].astype(str).to_numpy() if sample_key and sample_key in lineage.obs else np.repeat("single_object", lineage.n_obs)
    groups = lineage.obs[group_key].astype(str).to_numpy() if group_key and group_key in lineage.obs else np.repeat("all", lineage.n_obs)
    sample_summary = sample_cluster_table(samples, groups, labels, target, positive)
    highest_tests = highest_cluster_tests(sample_summary, highest)

    geometry = pd.DataFrame({
        "barcode": lineage.obs_names.astype(str), "UMAP1": coordinates[:, 0], "UMAP2": coordinates[:, 1],
        "cluster": labels, "sample": samples, "group": groups,
    })
    geometry_path = output / "02_shared_geometry" / f"{cfg['dataset_id']}_lineage_umap.tsv.gz"
    geometry.to_csv(geometry_path, sep="\t", index=False, compression="gzip")
    geometry_hash = hashlib.sha256(geometry[["barcode", "UMAP1", "UMAP2"]].to_csv(sep="\t", index=False).encode()).hexdigest()
    marker.to_csv(output / "05_tables" / "Figure2_marker_allowlist_summary.tsv", sep="\t", index=False)
    ledger.to_csv(output / "05_tables" / "Figure2_selected_marker_ledger.tsv", sep="\t", index=False)
    sample_summary.to_csv(output / "05_tables" / "Figure2_sample_cluster_summary.tsv", sep="\t", index=False)
    highest_tests.to_csv(output / "05_tables" / "Figure2_highest_target_cluster_tests_P_FDR.tsv", sep="\t", index=False)

    palette = dict(zip(cluster_order, plt.cm.tab20(np.linspace(0, 1, len(cluster_order)))))
    mode = cfg.get("panel_mode", "target_cluster")
    dpi = int(cfg.get("style", {}).get("dpi", 600))
    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.22, 1.0], left=0.075, right=0.935, bottom=0.06,
                           top=0.95, hspace=0.42, wspace=0.38)

    panel_a = fig.add_subplot(grid[0, 0])
    plot_cluster_umap(panel_a, coordinates, labels, palette, "A  Unsupervised lineage clusters")
    panel_a.legend(title="Cluster", loc="upper center", bbox_to_anchor=(0.5, -0.03), frameon=False,
                   markerscale=4, ncol=min(7, len(cluster_order)), columnspacing=0.7, handletextpad=0.4)

    if mode == "target_cluster":
        panel_b = fig.add_subplot(grid[0, 1])
        image = panel_b.scatter(coordinates[:, 0], coordinates[:, 1], s=1.2, c=target, cmap="viridis", linewidths=0,
                                rasterized=True)
        panel_b.set_title(f"B  {target_gene} expression", loc="left", fontweight="bold")
        clean_umap(panel_b)
        colorbar = fig.colorbar(image, ax=panel_b, fraction=0.045, pad=0.02)
        colorbar.set_label("log1p normalized expression")
        dot_axis = fig.add_subplot(grid[1, :])
        dot_letter = "C"
    else:
        dot_axis = fig.add_subplot(grid[0, 1])
        dot_letter = "B"

    dot = marker[marker.gene.isin(displayed_genes)]
    genes = displayed_genes[::-1]
    means = dot.pivot(index="gene", columns="cluster", values="mean_expression").reindex(index=genes, columns=cluster_order)
    percentages = dot.pivot(index="gene", columns="cluster", values="percent_expressed").reindex(index=genes, columns=cluster_order)
    zscore = means.sub(means.mean(axis=1), axis=0).div(means.std(axis=1).replace(0, np.nan), axis=0).fillna(0).clip(-2, 2)
    for y, gene in enumerate(genes):
        for x, cluster in enumerate(cluster_order):
            dot_axis.scatter(x, y, s=10 + 0.8 * percentages.loc[gene, cluster], c=zscore.loc[gene, cluster], cmap="RdBu_r",
                             vmin=-2, vmax=2, edgecolors="none")
    dot_axis.set_xticks(range(len(cluster_order)), cluster_order)
    dot_axis.set_yticks(range(len(genes)), genes)
    for tick in dot_axis.get_yticklabels():
        tick.set_fontstyle("italic")
        if tick.get_text() in mandatory:
            tick.set_color("#B6322C")
            tick.set_fontweight("bold")
    dot_axis.set_xlabel("Unsupervised cluster")
    dot_axis.set_title(f"{dot_letter}  Mandatory markers plus exact top 3 per cluster", loc="left", fontweight="bold")
    dot_axis.spines[["top", "right"]].set_visible(False)

    panel_contract = {"A": "unsupervised_lineage_umap", dot_letter: "fixed_allowlist_marker_dotplot"}
    if mode == "target_cluster":
        panel_d = fig.add_subplot(grid[2, 0])
        plot_cluster_umap(panel_d, coordinates, labels, palette, f"D  Highest {target_gene} cluster: {highest}", highlight=highest)
        panel_e = fig.add_subplot(grid[2, 1])
        co_summary = sample_summary.groupby("cluster").positive_fraction.agg(["mean", "sem"]).reindex(cluster_order)
        panel_e.bar(np.arange(len(cluster_order)), co_summary["mean"] * 100, color=[palette[x] for x in cluster_order],
                    edgecolor="#333333", linewidth=0.5)
        rng = np.random.default_rng(seed)
        for x, cluster in enumerate(cluster_order):
            values = sample_summary.loc[sample_summary.cluster.eq(cluster), "positive_fraction"].to_numpy() * 100
            panel_e.scatter(x + rng.uniform(-0.11, 0.11, len(values)), values, s=13, c="#555555", alpha=0.7, linewidths=0)
        panel_e.set_xticks(range(len(cluster_order)), cluster_order)
        panel_e.set_xlabel("Unsupervised cluster")
        panel_e.set_ylabel("+".join(coexpression) + "+ cells (%)")
        panel_e.set_title("E  Sample-unit coexpression summary", loc="left", fontweight="bold")
        panel_e.spines[["top", "right"]].set_visible(False)
        panel_contract.update({"B": "target_gene_feature_umap", "D": "highest_target_cluster_highlight", "E": "sample_unit_coexpression_by_cluster"})
    else:
        composition_axis = fig.add_subplot(grid[1, :])
        registry = sample_summary[["sample", "group"]].drop_duplicates()
        complete = registry.assign(_key=1).merge(
            pd.DataFrame({"cluster": cluster_order, "_key": 1}), on="_key").drop(columns="_key")
        composition = complete.merge(sample_summary[["sample", "group", "cluster", "n_cells"]],
                                     on=["sample", "group", "cluster"], how="left")
        composition["n_cells"] = composition["n_cells"].fillna(0).astype(int)
        totals = composition.groupby(["sample", "group"], as_index=False).n_cells.sum().rename(columns={"n_cells": "total"})
        composition = composition.merge(totals, on=["sample", "group"], how="left")
        composition["fraction"] = composition.n_cells / composition.total
        group_mean = composition.groupby(["group", "cluster"]).fraction.mean().unstack(fill_value=0).reindex(columns=cluster_order)
        bottom = np.zeros(len(group_mean))
        for cluster in cluster_order:
            composition_axis.bar(np.arange(len(group_mean)), group_mean[cluster] * 100, bottom=bottom * 100,
                                 color=palette[cluster], label=cluster, width=0.75)
            bottom += group_mean[cluster].to_numpy()
        composition_axis.set_xticks(range(len(group_mean)), group_mean.index, rotation=30, ha="right")
        composition_axis.set_ylabel("Mean sample-level composition (%)")
        composition_axis.set_title("C  Registered-group composition", loc="left", fontweight="bold")
        composition_axis.spines[["top", "right"]].set_visible(False)
        comparisons = cfg.get("comparisons", [])
        for index, letter in enumerate(["D", "E"]):
            axis = fig.add_subplot(grid[2, index])
            if index >= len(comparisons):
                axis.text(0.5, 0.5, "NOT_ESTIMABLE\nNo preregistered contrast", ha="center", va="center", fontweight="bold")
                axis.axis("off")
                panel_contract[letter] = "NOT_ESTIMABLE_missing_contrast"
                continue
            result = contrast_table(sample_summary, comparisons[index])
            result.to_csv(output / "05_tables" / f"Figure2{letter}_{comparisons[index].get('name', index + 1)}_P_FDR.tsv", sep="\t", index=False)
            y = np.arange(len(result))
            axis.axvline(0, color="#777777", linewidth=0.8)
            axis.barh(y, result.mean_fraction_difference * 100, color=[palette[x] for x in result.cluster], edgecolor="#333333", linewidth=0.5)
            axis.set_yticks(y, result.cluster)
            axis.set_xlabel("Mean proportion difference (percentage points)")
            axis.set_title(f"{letter}  {result.comparison.iloc[0]}", loc="left", fontweight="bold")
            for yy, row in enumerate(result.itertuples(index=False)):
                label = "NE" if not np.isfinite(row.Mann_Whitney_P) else f"P={row.Mann_Whitney_P:.3g}; FDR={row.BH_FDR:.3g}"
                axis.text(0.98, yy, label, transform=axis.get_yaxis_transform(), va="center", ha="right",
                          bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.5})
            axis.spines[["top", "right"]].set_visible(False)
            panel_contract[letter] = f"sample_level_group_contrast:{result.comparison.iloc[0]}"
        panel_contract["C"] = "registered_group_composition"

    fig.suptitle(str(cfg.get("display_title", f"{cfg['dataset_id']}: unsupervised subtype analysis")), x=0.075, y=0.987, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.018, "One frozen lineage UMAP is used throughout. Cells are not replicates; P and BH FDR use registered sample units where estimable.",
             ha="center", fontsize=max(8, float(cfg.get("style", {}).get("font_min_pt", 8))), fontstyle="italic")
    base = output / "04_figures" / f"Figure2_{cfg['dataset_id']}_A-E_{mode}_A4"
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".png"), dpi=dpi)
    plt.close(fig)

    manifest = {
        "dataset_id": cfg["dataset_id"], "input_h5ad": str(input_path), "input_sha256": sha256(input_path),
        "input_shape": list(source.shape), "lineage_shape": list(lineage.shape), "lineage_gate": cfg["lineage_gate"],
        "sample_key": sample_key, "group_key": group_key, "expression_source": expression_source,
        "embedding_source": representation, "seed": seed, "resolution": float(embedding_cfg.get("resolution", 0.55)),
        "n_clusters": len(cluster_order), "clusters": cluster_order, "target_gene": target_gene,
        "highest_target_cluster": highest, "coexpression_genes": coexpression, "positive_threshold": threshold,
        "panel_mode": mode, "panel_contract": panel_contract, "geometry_path": str(geometry_path),
        "geometry_sha256": geometry_hash, "old_global_X_umap_used": False,
        "canvas_mm": [210, 297], "font_family": cfg.get("style", {}).get("font_family", "Arial"),
        "font_range_pt": [cfg.get("style", {}).get("font_min_pt", 8), cfg.get("style", {}).get("font_max_pt", 12)],
    }
    (output / "00_manifest" / "Figure2_A-E_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    shutil.copy2(args.config, output / "00_manifest" / "figure2_config.json")
    shutil.copy2(Path(__file__).resolve(), output / "08_code" / Path(__file__).name)
    input_audit = {
        "input_h5ad": str(input_path), "sha256": sha256(input_path), "shape": list(source.shape),
        "obs_columns": list(map(str, source.obs.columns)), "layers": list(map(str, source.layers.keys())),
        "obsm_keys": list(map(str, source.obsm.keys())), "species": cfg.get("species"),
        "lineage_cells": int(lineage.n_obs), "lineage_genes": int(lineage.n_vars),
    }
    (output / "01_input_audit" / "input_audit.json").write_text(json.dumps(input_audit, ensure_ascii=False, indent=2), encoding="utf-8")
    audit = {
        "status": "PASS", "python": sys.version, "platform": platform.platform(), "anndata": ad.__version__, "scanpy": sc.__version__,
        "n_cells_input": source.n_obs, "n_cells_lineage": lineage.n_obs, "n_genes": lineage.n_vars,
        "mandatory_markers_present": all(gene in lineage.var_names for gene in mandatory),
        "exact_top3_per_cluster": bool((ledger.groupby("label").size() == 3).all()),
        "geometry_frozen": geometry_path.exists(), "old_global_X_umap_used": False,
        "p_and_fdr_table": str(output / "05_tables" / "Figure2_highest_target_cluster_tests_P_FDR.tsv"),
    }
    (output / "10_qa" / "Figure2_A-E_QA.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "07_reproduction" / "run_command.txt").write_text(f'python "{Path(__file__).resolve()}" --config "{args.config.resolve()}"\n', encoding="utf-8")
    relative_png = base.with_suffix(".png").relative_to(output).as_posix()
    relative_pdf = base.with_suffix(".pdf").relative_to(output).as_posix()
    entry = output / "06_html" / f"OPEN_ME_FIRST_{cfg['dataset_id']}_Figure2_complete.html"
    entry.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>Figure 2A-E</title>"
        "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:auto;padding:24px}"
        "h1{font-size:12pt}p,li{font-size:8pt;line-height:1.45}img{width:100%;height:auto;border:1px solid #ddd}</style>"
        f"</head><body><h1>{html.escape(str(cfg['dataset_id']))}: Figure 2A-E</h1>"
        f"<p>Mode: {html.escape(mode)}; frozen geometry SHA-256: <code>{geometry_hash}</code>.</p>"
        f"<p><a href='../{relative_pdf}'>Vector PDF</a> | "
        "<a href='../05_tables/Figure2_selected_marker_ledger.tsv'>marker ledger</a> | "
        "<a href='../10_qa/Figure2_A-E_QA.json'>QA</a> | "
        "<a href='../00_manifest/Figure2_A-E_manifest.json'>manifest</a></p>"
        f"<img src='../{relative_png}' alt='Figure 2A-E'></body></html>", encoding="utf-8"
    )
    output_rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file():
            output_rows.append({"relative_path": path.relative_to(output).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(output_rows).to_csv(output / "10_qa" / "output_index.tsv", sep="\t", index=False)
    print(json.dumps({"manifest": manifest, "qa": audit}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
