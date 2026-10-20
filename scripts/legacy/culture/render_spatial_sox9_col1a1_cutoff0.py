from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


OUT=Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
HUMAN_DATA=Path('data/external/INPUT_be81937a46/Wound_Healing_Annotation_Output_White')
HUMAN_BOUND=Path('data/external/INPUT_8665934387/mask')
RAT_DATA=Path('data/external/INPUT_01762fefba/spatial_h5ad_compatible')
RAT_BOUND=Path('data/external/INPUT_df0bc4e82e/annotation_output')
HUMAN=[("2 mph","2mph",HUMAN_DATA/"2mph_A03699G6.SCT_Annotated.h5ad"),("6 mph-1","6mph_p1",HUMAN_DATA/"6mph_A02988A2.SCT_part1_Annotated.h5ad"),("6 mph-2","6mph_p2",HUMAN_DATA/"6mph_A02988A2.SCT_part2_Annotated.h5ad"),("9 mph","9mph",HUMAN_DATA/"9mph_A03091C5.SCT_Annotated.h5ad")]
RAT=[("E16.5","6 h","E16.5_6H_C02938E1","Rat_E16.5_6H"),("E16.5","1 d","E16.5_1D_A03091D1","Rat_E16.5_1D"),("E16.5","2 d","E16.5_2D_A03091D4","Rat_E16.5_2D"),("E16.5","3 d","E16.5_3D_D03056D5","Rat_E16.5_3D"),("E16.5","6 d","E16.5_6D_C02938F4","Rat_E16.5_6D"),("Adult","1 d","Adult_1d_C01830C1D1","Rat_Adult_1d"),("Adult","2 d","Adult_2d_C01830E4F4","Rat_Adult_2d"),("Adult","5 d","Adult_5d1_C01830E6F6","Rat_Adult_5d"),("Adult","10 d","Adult_10d_A02988C5","Rat_Adult_10d"),("Adult","19 d","Adult_19d_B03424F4","Rat_Adult_19d")]


def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()


def dense(x): return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()
def rotate(xy,d):
    t=np.deg2rad(d); r=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    return np.asarray(xy,float)@r.T
def line_x(line,y):
    (x1,y1),(x2,y2)=np.asarray(line,float)
    return np.repeat((x1+x2)/2,len(y)) if np.isclose(y1,y2) else x1+(np.asarray(y)-y1)*(x2-x1)/(y2-y1)
def signed_distance(xy,line):
    p1,p2=np.asarray(line,float); v=p2-p1; rel=np.asarray(xy)-p1
    return (v[0]*rel[:,1]-v[1]*rel[:,0])/np.linalg.norm(v)
def signflip(d):
    d=np.asarray(d,float); d=d[np.isfinite(d)&(d!=0)]
    if not len(d): return np.nan
    obs=d.mean(); null=[]
    for mask in range(2**len(d)):
        s=np.array([1 if (mask>>i)&1 else -1 for i in range(len(d))]); null.append(np.mean(np.abs(d)*s))
    return float(np.mean(np.asarray(null)>=obs-1e-15))
def bh(p):
    p=np.asarray(p,float); out=np.full_like(p,np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    pv=p[ok]; o=np.argsort(pv); r=pv[o]; q=np.minimum.accumulate((r*len(r)/np.arange(1,len(r)+1))[::-1])[::-1]; inv=np.empty(len(o),int); inv[o]=np.arange(len(o)); out[ok]=np.clip(q[inv],0,1); return out
def style():
    mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":10,"axes.labelsize":9,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":.8})
def map_panel(ax,df,title,line):
    fib=df.is_fib.to_numpy(); pos=df.positive.to_numpy()
    # Exact old visual language: light-grey context, bright-red target, white canvas.
    ax.scatter(df.loc[fib,"x"],df.loc[fib,"y"],s=.35,c="#E0E0E0",alpha=.80,linewidths=0,rasterized=True)
    ax.scatter(df.loc[pos,"x"],df.loc[pos,"y"],s=1.5,c="#FF1A1A",alpha=.95,linewidths=0,rasterized=True)
    if line is not None:
        y0,y1=np.nanpercentile(df.loc[fib,"y"],[1,99]); ys=np.array([y0,y1]); ax.plot(line_x(line,ys),ys,"k--",lw=.9)
    ax.set_title(f"{title}  n={int(pos.sum())}",fontweight="bold",pad=2); ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")


def human():
    cells=[]; rows=[]; anchors={}; inputs=[]
    prepared=pd.read_csv(OUT/"03_analysis"/"human_spatial_SOX9_COL1A1_double_prepared.tsv.gz",sep="\t")
    prepared["is_fib"]=prepared.is_fib.astype(bool); prepared["positive"]=prepared.positive.astype(bool)
    anchor_cache=json.loads((OUT/"03_analysis"/"human_spatial_anchors_prepared.json").read_text(encoding="utf-8"))
    for display,key,path in HUMAN:
        ap=HUMAN_BOUND/f"{key}_spatial_anchors.json"; ainfo=anchor_cache[key]; anchors[key]=ainfo
        df=prepared[prepared.section==display].copy(); fib=df.is_fib.to_numpy(); pos=df.positive.to_numpy(); region=df.region.to_numpy(); cells.append(df)
        for reg in ["Inside scar","Outside scar"]:
            m=fib&(region==reg); rows.append({"section":display,"region":reg,"n_fib_bins":int(m.sum()),"n_positive":int(np.sum(m&pos)),"fraction":float(np.mean(pos[m])) if m.any() else np.nan})
        inputs.append({"section":display,"input":str(path),"sha256":sha256(path),"anchor":str(ap),"anchor_sha256":sha256(ap)})
    summary=pd.DataFrame(rows); wide=summary.pivot(index="section",columns="region",values="fraction").reindex([x[0] for x in HUMAN]); wide["inside_minus_outside"]=wide["Inside scar"]-wide["Outside scar"]
    p=signflip(wide.inside_minus_outside); fdr=p
    summary.to_csv(OUT/"05_tables"/"Figure1E_human_SOX9_COL1A1_gt0_by_section.tsv",sep="\t",index=False); wide.reset_index().to_csv(OUT/"05_tables"/"Figure1E_human_paired_effects_P_FDR.tsv",sep="\t",index=False)
    fig=plt.figure(figsize=(210/25.4,297/25.4)); gs=fig.add_gridspec(3,2,left=.07,right=.97,bottom=.07,top=.92,hspace=.28,wspace=.16,height_ratios=[1,1,.75])
    for i,((display,key,_),df) in enumerate(zip(HUMAN,cells)): map_panel(fig.add_subplot(gs[i//2,i%2]),df,display,anchors[key]["origin_line"])
    ax=fig.add_subplot(gs[2,0]);
    for sec,row in wide.iterrows(): ax.plot([0,1],[row["Outside scar"]*100,row["Inside scar"]*100],color="#777",lw=.8); ax.scatter([0,1],[row["Outside scar"]*100,row["Inside scar"]*100],c=["#BDBDBD","#FF1A1A"],edgecolor="black",lw=.4,s=25)
    ax.set_xticks([0,1],["Outside scar","Inside scar"]); ax.set_ylabel("SOX9+COL1A1+ Fib bins (%)"); ax.set_title("Section-level paired proportions",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False); ax.text(.02,.98,f"P={p:.4g}; BH FDR={fdr:.4g}",transform=ax.transAxes,va="top")
    ax=fig.add_subplot(gs[2,1]); d=wide.inside_minus_outside*100; ax.axhline(0,color="#777",lw=.8); ax.bar(range(len(d)),d,color="#FF1A1A",edgecolor="black",lw=.5); ax.set_xticks(range(len(d)),d.index); ax.set_ylabel("Inside - outside (percentage points)"); ax.set_title("Effect by section",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
    fig.suptitle(r"Human spatial SOX9$^{+}$COL1A1$^{+}$ fibroblasts",x=.07,y=.98,ha="left",fontsize=12,fontweight="bold"); fig.text(.07,.018,"Gate: source Fib annotation and normalized X >0 for SOX9 and COL1A1. Grey: other Fib bins; red: double-positive bins.\nStatistical unit: section/section-part.",fontsize=8)
    base=OUT/"04_figures"/"Figure1E_human_spatial_SOX9_COL1A1_gt0_oldstyle_A4"; fig.savefig(base.with_suffix(".pdf")); fig.savefig(base.with_suffix(".png"),dpi=600); plt.close(fig)
    rec={"gate":"Fib AND SOX9>0 AND COL1A1>0 in normalized source X","P_value":p,"BH_FDR":fdr,"inputs":inputs,"visual_style":{"context":"#E0E0E0","target":"#FF1A1A"},"canvas_mm":[210,297]}; (OUT/"00_manifest"/"Figure1E_human_spatial_manifest_v05.json").write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding="utf-8")


def rat():
    cells=[]; rows=[]; anchors={}; inputs=[]
    prepared=pd.read_csv(OUT/"03_analysis"/"rat_spatial_Sox9_Col1a1_double_prepared.tsv.gz",sep="\t")
    prepared["is_fib"]=prepared.is_fib.astype(bool); prepared["positive"]=prepared.positive.astype(bool)
    anchor_cache=json.loads((OUT/"03_analysis"/"rat_spatial_anchors_prepared.json").read_text(encoding="utf-8"))
    for stage,day,sample,akey in RAT:
        path=RAT_DATA/f"{sample}_Annotated_Spatial.h5ad"; ap=RAT_BOUND/f"{akey}_spatial_anchors.json"; info=anchor_cache[sample]; anchors[sample]=info
        df=prepared[prepared["sample"]==sample].copy(); fib=df.is_fib.to_numpy(); pos=df.positive.to_numpy(); zone=df.zone.to_numpy(); cells.append(df)
        for z in ["Wound-centre third","Flanking third"]:
            m=fib&(zone==z); rows.append({"stage":stage,"day":day,"sample":sample,"zone":z,"n_fib_bins":int(m.sum()),"n_positive":int(np.sum(m&pos)),"fraction":float(np.mean(pos[m])) if m.any() else np.nan})
        inputs.append({"sample":sample,"input":str(path),"sha256":sha256(path),"anchor":str(ap),"anchor_sha256":sha256(ap)})
    summary=pd.DataFrame(rows); wide=summary.pivot(index=["stage","day","sample"],columns="zone",values="fraction").reset_index(); wide["centre_minus_flank"]=wide["Wound-centre third"]-wide["Flanking third"]
    stat=[]
    for stage,q in wide.groupby("stage",sort=False): stat.append({"stage":stage,"n_sections":len(q),"mean_difference":q.centre_minus_flank.mean(),"one_sided_exact_signflip_P":signflip(q.centre_minus_flank)})
    stat=pd.DataFrame(stat); stat["BH_FDR_across_stages"]=bh(stat.one_sided_exact_signflip_P)
    summary.to_csv(OUT/"05_tables"/"Figure3G_rat_Sox9_Col1a1_gt0_by_section.tsv",sep="\t",index=False); wide.to_csv(OUT/"05_tables"/"Figure3G_rat_paired_effects.tsv",sep="\t",index=False); stat.to_csv(OUT/"05_tables"/"Figure3G_rat_stage_statistics_P_and_FDR.tsv",sep="\t",index=False)
    fig=plt.figure(figsize=(210/25.4,297/25.4)); gs=fig.add_gridspec(3,5,left=.045,right=.985,bottom=.08,top=.91,hspace=.28,wspace=.12,height_ratios=[1,1,.85])
    for i,((stage,day,sample,_),df) in enumerate(zip(RAT,cells)): map_panel(fig.add_subplot(gs[i//5,i%5]),df,f"{stage} {day}",anchors[sample]["origin_line"])
    for j,stage in enumerate(["E16.5","Adult"]):
        ax=fig.add_subplot(gs[2,j*2:(j+1)*2]); q=wide[wide.stage==stage]
        for _,row in q.iterrows():
            yy=[row["Flanking third"]*100,row["Wound-centre third"]*100]
            ax.plot([0,1],yy,color="#777",lw=.8)
            ax.scatter([0,1],yy,c=["#BDBDBD","#FF1A1A"],edgecolor="black",lw=.4,s=22)
        sr=stat[stat.stage==stage].iloc[0]; ax.set_xticks([0,1],["Flank","Centre"]); ax.set_ylabel("Sox9+Col1a1+ Fib bins (%)"); ax.set_title(f"{stage}: P={sr.one_sided_exact_signflip_P:.4g}; FDR={sr.BH_FDR_across_stages:.4g}",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
    ax=fig.add_subplot(gs[2,4]); d=wide.centre_minus_flank*100; ax.axhline(0,color="#777",lw=.8); ax.bar(range(len(d)),d,color=["#E28E2C" if x=="E16.5" else "#FF1A1A" for x in wide.stage],edgecolor="black",lw=.4); ax.set_xticks([]); ax.set_ylabel("Centre - flank (pp)"); ax.set_title("Section effects",loc="left",fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
    fig.suptitle(r"Rat spatial Sox9$^{+}$Col1a1$^{+}$ fibroblasts",x=.045,y=.98,ha="left",fontsize=12,fontweight="bold"); fig.text(.045,.018,"Gate: source Fib annotation and normalized X >0 for Sox9 and Col1a1. Grey: other Fib bins; red: double-positive bins.\nAnchors define wound centre versus flank, not histological scar boundaries.",fontsize=8)
    base=OUT/"04_figures"/"Figure3G_rat_spatial_Sox9_Col1a1_gt0_oldstyle_A4"; fig.savefig(base.with_suffix(".pdf")); fig.savefig(base.with_suffix(".png"),dpi=600); plt.close(fig)
    rec={"gate":"Fib AND Sox9>0 AND Col1a1>0 in normalized source X","boundary_semantics":"wound-centre versus flank; not validated scar inside/outside","statistics":stat.to_dict(orient="records"),"inputs":inputs,"visual_style":{"context":"#E0E0E0","target":"#FF1A1A"},"canvas_mm":[210,297]}; (OUT/"00_manifest"/"Figure3G_rat_spatial_manifest_v05.json").write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding="utf-8")


def main():
    for s in ["00_manifest","04_figures","05_tables"]:(OUT/s).mkdir(parents=True,exist_ok=True)
    style(); human(); rat(); print("spatial complete",flush=True)


if __name__=="__main__":main()
