from __future__ import annotations

import json
from pathlib import Path

import h5py
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


def decode_strings(values):
    return np.array([x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values], dtype=object)


def read_direct(path, genes):
    """Read only the required H5AD fields, bypassing unrelated newer /uns encodings."""
    with h5py.File(path, "r") as f:
        var_names = decode_strings(f["var"]["_index"][:])
        pos = {g: int(np.where(var_names == g)[0][0]) for g in genes}
        xg = f["X"]
        shape = tuple(int(x) for x in xg.attrs["shape"])
        x = sparse.csr_matrix((xg["data"][:], xg["indices"][:], xg["indptr"][:]), shape=shape)
        vals = {g: np.asarray(x[:, j].toarray()).ravel() for g, j in pos.items()}
        sg = f["obs"]["sub_labels"]
        categories = decode_strings(sg["categories"][:]); codes = sg["codes"][:].astype(int)
        sub = np.where(codes >= 0, categories[np.clip(codes, 0, len(categories)-1)], "NA")
        xy = np.asarray(f["obsm"]["spatial"][:], dtype=float)
    return sub, xy, vals
def rotate(xy,d):
    t=np.deg2rad(d); r=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    return np.asarray(xy,float)@r.T
def line_x(line,y):
    (x1,y1),(x2,y2)=np.asarray(line,float)
    return np.repeat((x1+x2)/2,len(y)) if np.isclose(y1,y2) else x1+(np.asarray(y)-y1)*(x2-x1)/(y2-y1)
def signed_distance(xy,line):
    p1,p2=np.asarray(line,float); v=p2-p1; rel=np.asarray(xy)-p1
    return (v[0]*rel[:,1]-v[1]*rel[:,0])/np.linalg.norm(v)


def main():
    d=OUT/"03_analysis"; d.mkdir(parents=True,exist_ok=True)
    human=[]; ha={}
    for display,key,path in HUMAN:
        info=json.loads((HUMAN_BOUND/f"{key}_spatial_anchors.json").read_text(encoding="utf-8")); ha[key]=info
        sub,xy,vals=read_direct(path,["SOX9","COL1A1"]); fib=np.char.startswith(sub.astype(str),"Fib_"); xy=rotate(xy,float(info.get("rotation_applied",0)))
        pos=fib.copy()
        for g in vals: pos &= vals[g]>0.0
        region=np.where(xy[:,0]>line_x(info["origin_line"],xy[:,1]),"Outside scar","Inside scar")
        human.append(pd.DataFrame({"section":display,"key":key,"x":xy[:,0],"y":xy[:,1],"is_fib":fib,
                                   "SOX9":vals["SOX9"],"COL1A1":vals["COL1A1"],"positive":pos,"region":region}))
        print("human",display,flush=True)
    pd.concat(human,ignore_index=True).to_csv(d/"human_spatial_SOX9_COL1A1_double_prepared.tsv.gz",sep="\t",index=False,compression="gzip")
    (d/"human_spatial_anchors_prepared.json").write_text(json.dumps(ha,ensure_ascii=False,indent=2),encoding="utf-8")
    rat=[]; ra={}
    for stage,day,sample,akey in RAT:
        info=json.loads((RAT_BOUND/f"{akey}_spatial_anchors.json").read_text(encoding="utf-8")); ra[sample]=info; path=RAT_DATA/f"{sample}_Annotated_Spatial.h5ad"
        sub,xy,vals=read_direct(path,["Sox9","Col1a1"]); fib=np.char.startswith(sub.astype(str),"Fib_"); xy=rotate(xy,float(info.get("rotation_applied",0)))
        pos=fib.copy()
        for g in vals: pos &= vals[g]>0.0
        dist=np.abs(signed_distance(xy,info["origin_line"])); q1,q2=np.quantile(dist[fib],[1/3,2/3]); zone=np.full(len(fib),"Middle third",object); zone[dist<=q1]="Wound-centre third"; zone[dist>=q2]="Flanking third"
        rat.append(pd.DataFrame({"stage":stage,"day":day,"sample":sample,"x":xy[:,0],"y":xy[:,1],"is_fib":fib,
                                 "Sox9":vals["Sox9"],"Col1a1":vals["Col1a1"],"positive":pos,"zone":zone}))
        print("rat",sample,flush=True)
    pd.concat(rat,ignore_index=True).to_csv(d/"rat_spatial_Sox9_Col1a1_double_prepared.tsv.gz",sep="\t",index=False,compression="gzip")
    (d/"rat_spatial_anchors_prepared.json").write_text(json.dumps(ra,ensure_ascii=False,indent=2),encoding="utf-8")


if __name__=="__main__":main()
