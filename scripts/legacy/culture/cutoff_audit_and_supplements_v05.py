from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

OUT = Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
DATASETS = {
    "human_scRNA": {
        "h5ad": Path('data/external/INPUT_decd995edd/pbmc_final.h5ad'),
        "geometry": OUT / "03_analysis/human_Figure2_AE_composition/02_shared_geometry/human_Fib_cutoff0_lineage_umap.tsv.gz",
        "genes": ("SOX9", "COL1A1"), "normal": "Normal", "title": "Human Fib-lineage",
    },
    "rat_scRNA": {
        "h5ad": Path('data/external/INPUT_a5f4df094a/pbmc_final.raw_counts.h5ad'),
        "geometry": OUT / "03_analysis/rat_Figure3G_AE_composition/02_shared_geometry/rat_Fib_cutoff0_lineage_umap.tsv.gz",
        "genes": ("Sox9", "Col1a1"), "normal": "Normal", "title": "Rat Fib-lineage",
    },
}

def bh(values):
    p = np.asarray(values, float); out = np.full(p.shape, np.nan); ok = np.isfinite(p)
    if not ok.any(): return out
    pv = p[ok]; order = np.argsort(pv); ranked = pv[order]
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty(len(pv)); restored[order] = np.minimum(q, 1); out[ok] = restored
    return out

def normalized_two_genes(path: Path, barcodes, genes, target_sum=10000.0):
    a = ad.read_h5ad(path, backed="r")
    rows = pd.Index(a.obs_names.astype(str)).get_indexer(pd.Index(barcodes.astype(str)))
    if np.any(rows < 0): raise ValueError(f"Missing {int(np.sum(rows < 0))} barcodes in {path}")
    positions = pd.Index(a.var_names.astype(str)).get_indexer(genes)
    if np.any(positions < 0): raise ValueError(f"Missing genes {genes} in {path}")
    values = np.zeros((len(rows), 2), float)
    for start in range(0, len(rows), 2048):
        stop = min(start + 2048, len(rows)); rr = rows[start:stop]
        block = a.X[rr, :]
        if not sparse.issparse(block): block = sparse.csr_matrix(block)
        totals = np.asarray(block.sum(axis=1)).ravel()
        raw = np.asarray(block[:, positions].toarray(), float)
        scale = np.divide(target_sum, totals, out=np.zeros_like(totals, float), where=totals > 0)
        values[start:stop] = np.log1p(raw * scale[:, None])
    a.file.close()
    return values

def setup_style():
    mpl.rcParams.update({"font.family":"Arial", "font.size":8, "axes.titlesize":10, "axes.labelsize":9,
                         "xtick.labelsize":8, "ytick.labelsize":8, "legend.fontsize":8,
                         "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none", "axes.linewidth":.8})

def summarize(frame, dataset):
    rows=[]
    for dims in [["sample"], ["group"], ["cluster"], ["sample","group","cluster"]]:
        for keys, q in frame.groupby(dims, dropna=False, sort=False):
            if not isinstance(keys, tuple): keys=(keys,)
            row={"dataset":dataset,"level":"+".join(dims),"n_fib_cells":len(q),
                 "n_gene1_gt0":int((q.gene1>0).sum()),"n_gene2_gt0":int((q.gene2>0).sum()),
                 "n_double_gt1":int(q.double_gt1.sum()),"n_double_gt0":int(q.double_gt0.sum()),
                 "fraction_double_gt1":float(q.double_gt1.mean()),"fraction_double_gt0":float(q.double_gt0.mean())}
            row.update(dict(zip(dims,keys))); rows.append(row)
    return pd.DataFrame(rows)

def supplemental_plot(frame, info, highest, dataset):
    genes=info["genes"]; clusters=sorted(frame.cluster.unique(),key=lambda x:int(x[1:])); palette=dict(zip(clusters,plt.cm.tab20(np.linspace(0,1,len(clusters)))))
    fig=plt.figure(figsize=(210/25.4,297/25.4)); gs=fig.add_gridspec(2,2,left=.08,right=.96,bottom=.12,top=.91,hspace=.30,wspace=.28)
    ax=fig.add_subplot(gs[0,0]); im=ax.scatter(frame.UMAP1,frame.UMAP2,c=frame.gene1,s=1.2,cmap="viridis",linewidths=0,rasterized=True); ax.set_title(f"A  {genes[0]} expression",loc="left",fontweight="bold"); ax.axis("off"); cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.02); cb.set_label("log1p normalized expression")
    ax=fig.add_subplot(gs[0,1]); bg=frame.cluster.ne(highest); ax.scatter(frame.loc[bg,"UMAP1"],frame.loc[bg,"UMAP2"],s=1.2,c="#E4E4E4",alpha=.65,linewidths=0,rasterized=True); q=frame[~bg]; ax.scatter(q.UMAP1,q.UMAP2,s=2,c=palette[highest],linewidths=0,rasterized=True); ax.set_title(f"B  Highest observed {genes[0]} cluster: {highest}",loc="left",fontweight="bold"); ax.axis("off")
    sample=frame.groupby(["sample","group","cluster"],as_index=False).double_gt0.mean(); ax=fig.add_subplot(gs[1,0]); mean=sample.groupby("cluster").double_gt0.mean().reindex(clusters)*100; ax.bar(range(len(clusters)),mean,color=[palette[x] for x in clusters],edgecolor="#333",lw=.5); rng=np.random.default_rng(20260917)
    for x,c in enumerate(clusters):
        vals=sample.loc[sample.cluster.eq(c),"double_gt0"].to_numpy()*100; ax.scatter(x+rng.uniform(-.11,.11,len(vals)),vals,s=13,c="#555",alpha=.7,linewidths=0)
    ax.set_xticks(range(len(clusters)),clusters,rotation=45,ha="right"); ax.set_ylabel(f"{genes[0]}+{genes[1]}+ cells (%)"); ax.set_title("C  Sample-level cutoff 0 coexpression",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
    ax=fig.add_subplot(gs[1,1]); focus=sample[sample.cluster.eq(highest)].copy(); groups=list(pd.unique(focus.group)); groups=sorted(groups,key=lambda x:(x!=info["normal"],x)); rng=np.random.default_rng(20260917); data=[focus.loc[focus.group.eq(g),"double_gt0"].to_numpy()*100 for g in groups]; ax.boxplot(data,positions=range(len(groups)),widths=.55,showfliers=False,patch_artist=True,boxprops={"facecolor":"#F2F2F2"},medianprops={"color":"#B22222"})
    for x,vals in enumerate(data): ax.scatter(x+rng.uniform(-.08,.08,len(vals)),vals,s=18,c="#555",alpha=.75,linewidths=0)
    rows=[]; normal=info["normal"]
    for g in groups:
        if g==normal: continue
        a=np.asarray(data[groups.index(g)]); b=np.asarray(data[groups.index(normal)]); p=float(stats.mannwhitneyu(a,b,alternative="two-sided").pvalue) if len(a)>=3 and len(b)>=3 else np.nan
        rows.append({"dataset":dataset,"cluster":highest,"comparison":f"{g}_vs_{normal}","n_numerator":len(a),"n_denominator":len(b),"mean_difference_percentage_points":float(np.mean(a)-np.mean(b)),"Mann_Whitney_P":p})
    tests=pd.DataFrame(rows); tests["BH_FDR"]=bh(tests.Mann_Whitney_P) if len(tests) else []
    ax.set_xticks(range(len(groups)),groups,rotation=30,ha="right"); ax.set_ylabel(f"{genes[0]}+{genes[1]}+ cells in {highest} (%)"); ax.set_title(f"D  {highest} sample-level proportions",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
    if len(tests): ax.text(.02,.98,"; ".join(f"{r.comparison}: P={r.Mann_Whitney_P:.3g}, FDR={r.BH_FDR:.3g}" for r in tests.itertuples()),transform=ax.transAxes,va="top",fontsize=8,wrap=True)
    fig.suptitle(f"{info['title']}: {genes[0]}+{genes[1]}+ cutoff 0 supplements",x=.08,y=.98,ha="left",fontsize=12,fontweight="bold"); fig.text(.08,.035,"These are supplementary panels and do not replace the formal composition_contrasts Figure A-E. Cells are not inferential replicates.",fontsize=8,fontstyle="italic")
    base=OUT/"04_figures"/f"{dataset}_target_highlight_and_cutoff0_supplement_A4"; fig.savefig(base.with_suffix(".pdf")); fig.savefig(base.with_suffix(".svg")); fig.savefig(base.with_suffix(".png"),dpi=600); plt.close(fig)
    tests.to_csv(OUT/"05_tables"/f"{dataset}_highest_cluster_cutoff0_group_tests_P_FDR.tsv",sep="\t",index=False)

def scrna():
    totals=[]; all_summaries=[]
    for dataset,info in DATASETS.items():
        geom=pd.read_csv(info["geometry"],sep="\t"); vals=normalized_two_genes(info["h5ad"],geom.barcode,info["genes"])
        frame=geom.copy(); frame["gene1"]=vals[:,0]; frame["gene2"]=vals[:,1]; frame["double_gt0"]=(vals[:,0]>0)&(vals[:,1]>0); frame["double_gt1"]=(vals[:,0]>1)&(vals[:,1]>1)
        added=frame.double_gt0 & ~frame.double_gt1
        totals.append({"dataset":dataset,"expression_source":"X->normalize_total(10000)->log1p","gene1":info["genes"][0],"gene2":info["genes"][1],"n_fib":len(frame),"n_gene1_gt0":int((frame.gene1>0).sum()),"n_gene2_gt0":int((frame.gene2>0).sum()),"n_double_gt1":int(frame.double_gt1.sum()),"n_double_gt0":int(frame.double_gt0.sum()),"n_added":int(added.sum()),"added_gene1_min":float(frame.loc[added,"gene1"].min()) if added.any() else np.nan,"added_gene1_max":float(frame.loc[added,"gene1"].max()) if added.any() else np.nan,"added_gene2_min":float(frame.loc[added,"gene2"].min()) if added.any() else np.nan,"added_gene2_max":float(frame.loc[added,"gene2"].max()) if added.any() else np.nan})
        summary=summarize(frame,dataset); all_summaries.append(summary); summary.to_csv(OUT/"05_tables"/f"{dataset}_cutoff1_vs_cutoff0_by_sample_group_cluster.tsv",sep="\t",index=False)
        frame[["barcode","cluster","sample","group","gene1","gene2","double_gt1","double_gt0"]].to_csv(OUT/"05_tables"/f"{dataset}_cutoff_cell_audit.tsv.gz",sep="\t",index=False,compression="gzip")
        highest=frame.groupby("cluster").gene1.mean().sort_values(ascending=False).index[0]; supplemental_plot(frame,info,highest,dataset)
    pd.DataFrame(totals).to_csv(OUT/"05_tables"/"scRNA_cutoff1_vs_cutoff0_total_audit.tsv",sep="\t",index=False)
    pd.concat(all_summaries,ignore_index=True).to_csv(OUT/"05_tables"/"scRNA_cutoff1_vs_cutoff0_all_levels.tsv",sep="\t",index=False)
    return totals

def spatial():
    totals=[]
    h=pd.read_csv(OUT/"03_analysis/human_spatial_SOX9_COL1A1_double_prepared.tsv.gz",sep="\t"); h=h[h.is_fib.astype(bool)].copy(); h["double_gt0"]=(h.SOX9>0)&(h.COL1A1>0); h["double_gt1"]=(h.SOX9>1)&(h.COL1A1>1)
    for sec,q in h.groupby("section",sort=False):
        added=q.double_gt0&~q.double_gt1; totals.append({"dataset":"human_spatial","unit":sec,"n_fib":len(q),"n_double_gt1":int(q.double_gt1.sum()),"n_double_gt0":int(q.double_gt0.sum()),"n_added":int(added.sum()),"added_gene1_min":float(q.loc[added,"SOX9"].min()) if added.any() else np.nan,"added_gene1_max":float(q.loc[added,"SOX9"].max()) if added.any() else np.nan,"added_gene2_min":float(q.loc[added,"COL1A1"].min()) if added.any() else np.nan,"added_gene2_max":float(q.loc[added,"COL1A1"].max()) if added.any() else np.nan})
    hs=h.groupby(["section","region"],as_index=False).agg(n_fib=("double_gt0","size"),n_double_gt1=("double_gt1","sum"),n_double_gt0=("double_gt0","sum")); hs.to_csv(OUT/"05_tables"/"human_spatial_cutoff1_vs_cutoff0_by_section_region.tsv",sep="\t",index=False)
    r=pd.read_csv(OUT/"03_analysis/rat_spatial_Sox9_Col1a1_double_prepared.tsv.gz",sep="\t"); r=r[r.is_fib.astype(bool)].copy(); r["double_gt0"]=(r.Sox9>0)&(r.Col1a1>0); r["double_gt1"]=(r.Sox9>1)&(r.Col1a1>1)
    for sample,q in r.groupby("sample",sort=False):
        added=q.double_gt0&~q.double_gt1; totals.append({"dataset":"rat_spatial","unit":sample,"n_fib":len(q),"n_double_gt1":int(q.double_gt1.sum()),"n_double_gt0":int(q.double_gt0.sum()),"n_added":int(added.sum()),"added_gene1_min":float(q.loc[added,"Sox9"].min()) if added.any() else np.nan,"added_gene1_max":float(q.loc[added,"Sox9"].max()) if added.any() else np.nan,"added_gene2_min":float(q.loc[added,"Col1a1"].min()) if added.any() else np.nan,"added_gene2_max":float(q.loc[added,"Col1a1"].max()) if added.any() else np.nan})
    rs=r.groupby(["stage","sample","zone"],as_index=False).agg(n_fib=("double_gt0","size"),n_double_gt1=("double_gt1","sum"),n_double_gt0=("double_gt0","sum")); rs.to_csv(OUT/"05_tables"/"rat_spatial_cutoff1_vs_cutoff0_by_section_zone.tsv",sep="\t",index=False)
    pd.DataFrame(totals).to_csv(OUT/"05_tables"/"spatial_cutoff1_vs_cutoff0_total_audit.tsv",sep="\t",index=False)
    return totals

def main():
    setup_style(); sc=scrna(); sp=spatial(); audit={"status":"PASS","gate_human":"SOX9>0 & COL1A1>0","gate_rat":"Sox9>0 & Col1a1>0","old_cutoff":1.0,"new_cutoff":0.0,"scrna":sc,"spatial":sp}
    (OUT/"10_qa"/"CUTOFF_CHANGE_AUDIT.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8"); print(json.dumps(audit,ensure_ascii=False,indent=2))

if __name__=="__main__": main()
