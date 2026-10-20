from pathlib import Path
import json, sys, hashlib
import numpy as np, pandas as pd, openpyxl
from prepare_full_selected import ROOT, PROJECT, RNA
sys.stdout.reconfigure(encoding='utf8')
source=PROJECT/'20261009_TableS1_Figure1_4_actual_inputs_v05/outputs/Table_S1_Figure1_4_actual_scRNA_samples_v05.xlsx'
w=openpyxl.load_workbook(source,read_only=True,data_only=True)
it=w['Table S1'].iter_rows(values_only=True);columns=list(next(it));old=pd.DataFrame(list(it),columns=columns);w.close()
cells=pd.read_csv(ROOT/'tables/Figure1_selected19_full_QC_skin_cell_metadata.tsv.gz',sep='\t',index_col=0,low_memory=False,float_precision='round_trip')
fib=pd.read_csv(ROOT/'tables/Selected19_full_Fib_library_expression.csv',float_precision='round_trip').set_index('record_id')
public=pd.read_csv(PROJECT/'20261009_Figure1BC_selected7_controls_v02/audit/selected_cells_metadata_expression.tsv.gz',sep='\t',index_col=0,low_memory=False)
public=public.loc[~public.organ.eq('skin')]
ids=set(public.record_id)|set(cells.record_id)|set(old.loc[old['Record ID'].str.startswith('HFB'),'Record ID'])
s1=old.loc[old['Record ID'].isin(ids)].copy();assert len(s1)==104 and s1['Record ID'].is_unique
assert len(set(public.record_id))==77 and len(set(cells.record_id))==19
counts=cells.groupby('record_id').size()
for rid in fib.index:
    idx=s1.index[s1['Record ID'].eq(rid)][0];f=fib.loc[rid]
    s1.loc[idx,'Used panels']='Figure1B; Figure1C; Figure2A/B'
    s1.loc[idx,'Dataset']='Self-cohort human skin, selected full libraries'
    s1.loc[idx,'Figure1B cells']=int(counts[rid]);s1.loc[idx,'Figure1C target']='All Fib (original Fibroblast annotation)'
    s1.loc[idx,'Figure1C target cells']=int(f.Fib_cells)
    s1.loc[idx,'Figure1C SOX9 expression sum']=float(f.SOX9_mean*f.Fib_cells)
    s1.loc[idx,'Figure1C SOX9 mean expression']=float(f.SOX9_mean)
    s1.loc[idx,'Figure1C group']=f.condition_pair;s1.loc[idx,'Figure2A/B Fib cells']=int(f.Fib_cells)
    s1.loc[idx,'Figure2 clinical group']=f.clinical_group;s1.loc[idx,'Patient ID as recorded']=f.patient_id_source
    s1.loc[idx,'Source panel']='Selected full human skin cohort, 19 libraries'
    s1.loc[idx,'Biological replicate status']='Patient IDs retained from approved clinical mapping; repeated libraries pooled by patient for exploratory statistics'
    warning='User-selected after SOX9 expression review; exploratory cohort. Full original Fib gate, not reference-sampled Fib47.'
    if f.Fib_cells<20:warning+=' Fewer than 20 Fib cells; expression estimate is imprecise.'
    s1.loc[idx,'Warnings']=warning
for rid in public.record_id.unique():
    x=s1.loc[s1['Record ID'].eq(rid)].iloc[0]
    assert int(x['Figure1B cells'])==int(public.record_id.eq(rid).sum())
assert s1.loc[s1['Figure1B cells'].gt(0),'Figure1B cells'].sum()==465450
assert s1['Figure2A/B Fib cells'].sum()==20270 and s1['Figure4 all Fib cells'].sum()==59854
assert not s1['True tissue'].astype(str).str.fullmatch('PBMC|BM|Bone marrow',case=False).any()
s1=s1.sort_values('Record ID').reset_index(drop=True)
s1.to_csv(ROOT/'tables/Table_S1_current_selected19_Figure1_4_actual_scRNA.csv',index=False,encoding='utf-8-sig')
defs=pd.DataFrame([
    ['Scope','104 actual scRNA library records: 96 used by current Figure1B-C (19 skin + 77 other-organ records), plus 8 unchanged Figure4 cultured Fib libraries. Figure2 uses the same 19 skin library IDs. Spatial sections are not counted as scRNA libraries.'],
    ['Skin cells','208,532 full QC-pass skin cells; 20,270 original Fibroblast cells. No reference sampling or gene-positive prefilter.'],
    ['Skin target','All Fib, defined by original bulk_labels=Fibroblast. Fib_Myo is not the main skin target.'],
    ['Expression','Mean existing nonnegative source normalized RNA X within all gated Fib, including zeros. Public source expressions and target gates remain unchanged.'],
    ['QC','Same previous thresholds: source analysis_qc_pass, 200–8000 genes, 500–60000 counts, mt<=20%, HB<=1%.'],
    ['Statistical units','Skin: equal-library-weighted means within recorded patient IDs, 5 normal and 10 scar patient IDs. These IDs are retained from the clinical map, not independently revalidated identities.'],
    ['Selection limitation','The 19-library cohort was chosen after viewing SOX9 expression summaries. Its P/FDR results are exploratory and do not constitute independently prespecified confirmation.'],
    ['Figure1 geometry','Fresh joint PCA/neighbors/UMAP on all 465,450 selected/current cells. Human/mouse symbol-harmonized cross-organ visualization is descriptive.'],
    ['Figure2 geometry','Fresh full-lineage HVG/PCA/Leiden/UMAP from the immutable universal skill executable. C labels are unsupervised; highest SOX9 mean is an observed ranking, not a definitive subtype assignment.'],
    ['Culture library boundary','All eight Figure4 libraries are cultured skin-derived Fib treatment libraries. Donor and independent culture replicate IDs are not recorded. Two library labels per arm do not establish two biological replicates.'],
    ['No absent target as zero','SOX9 zeros within existing Fib are retained. Missing target populations elsewhere remain unavailable; no biological expression is imputed.'],
],columns=['Definition','Details'])
payload=[]
for name,d in [('Table S1',s1),('Definitions',defs)]:payload.append({'name':name,'headers':list(d.columns),'rows':json.loads(d.replace([np.inf,-np.inf],np.nan).to_json(orient='values',force_ascii=False))})
(ROOT/'audit/TableS1_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'audit/TableS1_scope.json').write_text(json.dumps({'source_workbook':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'current_scRNA_records':104,'Figure1_records':96,'Figure1_skin_records':19,'Figure1_non_skin_records':77,'Figure2_records':19,'Figure4_records':8,'PBMC_specimen_records':0,'Figure1_cells':465450,'Figure2_Fib_cells':20270,'Figure4_cells':59854,'source_workbook_not_overwritten':True},ensure_ascii=False,indent=2),encoding='utf8')
print('Table S1 payload',len(s1),'rows; no PBMC specimens',flush=True)
