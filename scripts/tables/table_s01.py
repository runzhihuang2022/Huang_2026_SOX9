#!/usr/bin/env python3
"""Export the already-reconciled, privacy-minimized actual-library Table S1."""
from pathlib import Path
import argparse
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[2];out=Path(a.out)
    if out.exists():raise FileExistsError('Refusing to overwrite '+str(out))
    d=pd.read_csv(root/'data/metadata/Table_S1_current_public_summary.csv',float_precision='round_trip')
    used=d['Figure1B cells'].gt(0)|d['Figure2A/B Fib cells'].gt(0)|d['Figure4 all Fib cells'].gt(0)
    assert used.all() and len(d)==104 and d['Record ID'].is_unique
    assert not d['True tissue'].str.contains('PBMC',case=False).any()
    assert d['Figure1B cells'].sum()==465450 and d['Figure2A/B Fib cells'].sum()==20270
    out.parent.mkdir(parents=True,exist_ok=True);d.to_csv(out,index=False,encoding='utf-8-sig')

if __name__=='__main__':main()
