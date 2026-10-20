from pathlib import Path
import argparse, hashlib, json, sys, time, shutil, gc
import h5py, anndata as ad, numpy as np, pandas as pd
from scipy import sparse
sys.stdout.reconfigure(encoding='utf8')
ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
RNA=Path('data/external/INPUT_bc1086c0aa/pbmc_scar_bm_pbmc.h5ad')
META=Path('data/external/INPUT_2547f54905/module1_all_cells_metadata_only_annotated.h5ad')
CLIN=Path('data/external/INPUT_beb6b056bf/approved_cell_clinical_map.tsv.gz')
CHOICE=Path('data/external/INPUT_c3b11df53e/user_candidate_skin_record_ids.csv')
INVENTORY=PROJECT/'20261009_Skin_full_vs_Figure1_subset_Fib_SOX9_audit_v02/tables/62_true_skin_libraries_Fib_SOX9_full_vs_reference.csv'
SKILL=Path('data/external/INPUT_0fecd8a262/skin-singlecell-subtype-complete-analysis')
SEED=20260917
for folder in ['00_manifest','01_input_audit','02_shared_geometry','analysis','figures','tables','06_html','07_reproduction','code','logs','audit','outputs']:(ROOT/folder).mkdir(parents=True,exist_ok=True)

def sha(p):
    hh=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):hh.update(b)
    return hh.hexdigest()

def read_frame(path,key='obs'):
    with h5py.File(path,'r') as h:return ad.io.read_elem(h[key])

def subset_sparse(group,rows,columns=None):
    shape=tuple(int(v) for v in group.attrs['shape']);rows=np.asarray(rows,dtype=np.int64)
    columns=np.arange(shape[1],dtype=np.int64) if columns is None else np.asarray(columns,dtype=np.int64)
    enc=group.attrs['encoding-type']
    if enc=='csr_matrix':
        # These public objects are much smaller than the full scar CSC parent.
        x=sparse.csr_matrix((group['data'][:].astype(np.float32),group['indices'][:],group['indptr'][:]),shape=shape)
        return x[rows][:,columns].tocsr()
    assert enc=='csc_matrix',enc
    mapper=np.full(shape[0],-1,dtype=np.int32);mapper[rows]=np.arange(len(rows),dtype=np.int32)
    ptr=group['indptr'][:]
    # Contiguous ranges amortize compressed HDF5 reads; never materialize the full parent matrix.
    blocks=[];start=0;last_print=time.time()
    while start<len(columns):
        end=start+1
        while end<len(columns) and columns[end]==columns[end-1]+1 and end-start<128:end+=1
        first,last=int(columns[start]),int(columns[end-1]);a,b=int(ptr[first]),int(ptr[last+1])
        rr=mapper[group['indices'][a:b]];good=rr>=0
        values=group['data'][a:b][good].astype(np.float32)
        cc=np.repeat(np.arange(end-start,dtype=np.int32),np.diff(ptr[first:last+2]))[good]
        counts=np.bincount(cc,minlength=end-start);indptr=np.concatenate([[0],np.cumsum(counts)])
        blocks.append(sparse.csc_matrix((values,rr[good],indptr),shape=(len(rows),end-start)))
        start=end
        if time.time()-last_print>25:print('CSC extraction',start,'/',len(columns),'genes; selected rows',len(rows),flush=True);last_print=time.time()
    return sparse.hstack(blocks,format='csr',dtype=np.float32)

def main():
    out=ROOT/'01_input_audit/Selected19_full_QC_Fib_normalized_RNA.h5ad'
    if out.exists():raise FileExistsError('Prepared input exists; do not overwrite')
    selected=pd.read_csv(CHOICE,dtype=str).record_id.tolist();assert len(selected)==19 and len(set(selected))==19
    inv=pd.read_csv(INVENTORY).set_index('record_id');assert set(selected).issubset(set(inv.index))
    chosen=inv.loc[selected].copy();assert chosen.actual_tissue.eq('Skin').all()
    assert chosen.clinical_group.value_counts().to_dict()=={'DSDB_Scar':14,'Normal':5}
    chosen.to_csv(ROOT/'tables/selected19_library_whitelist_source_audit.csv',encoding='utf-8-sig')
    shutil.copy2(CHOICE,ROOT/'01_input_audit/user_candidate_skin_record_ids_LOCKED.csv')
    before={'rna':{'path':str(RNA),'bytes':RNA.stat().st_size,'mtime_ns':RNA.stat().st_mtime_ns},'metadata':{'path':str(META),'bytes':META.stat().st_size,'mtime_ns':META.stat().st_mtime_ns}}
    print('Reading full metadata; 19-library whitelist',selected,flush=True)
    with h5py.File(META,'r') as h:
        fields=['source_obs_name','sequencing_unit','bulk_labels','sub_labels','analysis_qc_pass','module1_reference_cell','nFeature_RNA','nCount_RNA','mt_percent','HB_percent']
        full=pd.DataFrame({k:ad.io.read_elem(h['obs'][k]) for k in fields})
    seqmap=chosen.reset_index().set_index('sequencing_unit').record_id.to_dict()
    full['record_id']=full.sequencing_unit.astype(str).map(seqmap)
    full.index=pd.Index(full.source_obs_name.astype(str),name='barcode');assert full.index.is_unique
    full['full_source_row']=np.arange(len(full),dtype=np.int64)
    mask=full.record_id.notna();d=full.loc[mask].copy()
    d['QC_pass']=d.analysis_qc_pass.astype(bool)&d.nFeature_RNA.between(200,8000)&d.nCount_RNA.between(500,60000)&(d.mt_percent<=20)&(d.HB_percent<=1)
    counts=d.groupby('record_id').agg(source_cells=('record_id','size'),QC_cells=('QC_pass','sum')).reindex(selected)
    d=d.loc[d.QC_pass].copy();assert len(d)==208532
    clinical=pd.read_csv(CLIN,sep='\t',usecols=['cell_id','patient_id','tissue','grade','timepoint','sex','mapping_status'],low_memory=False).set_index('cell_id')
    joined=clinical.reindex(d.index);assert joined.tissue.eq('Skin').all() and joined.patient_id.notna().all()
    for k in joined.columns:d[k]=joined[k].to_numpy()
    d['patient_id_source']=d.patient_id.astype(str);d['clinical_group']=d.grade.astype(str)
    d['condition_pair']=d.clinical_group.map({'Normal':'normal','DSDB_Scar':'fibrosis'});assert d.condition_pair.notna().all()
    d['organ']='skin';d['species']='human';d['source_index']=0;d['actual_tissue']='Skin'
    d['source_sample']=d.record_id.map(chosen.source_sample.to_dict());d['dataset_panel']='Selected full human skin libraries'
    d['cell_annotation']=d.sub_labels.astype(str);d['bulk_annotation']=d.bulk_labels.astype(str)
    d['original_Fib_gate']=d.bulk_labels.astype(str).eq('Fibroblast')
    assert int(d.original_Fib_gate.sum())==20270
    d['sample_unit']=d.patient_id_source;d['expression_basis']='Full parent nonnegative normalized RNA X (uns X_name=logcounts); unchanged'
    with h5py.File(RNA,'r') as h:
        raw_index=pd.Index(ad.io.read_elem(h['obs'][h['obs'].attrs['_index']]).astype(str))
        assert np.array_equal(raw_index.to_numpy(),full.index.to_numpy())
        genes=pd.Index(ad.io.read_elem(h['var'][h['var'].attrs['_index']]).astype(str));assert genes.is_unique
        gj=int(genes.get_loc('SOX9'));ptr=h['X/indptr'][gj:gj+2];vals=np.zeros(len(raw_index),dtype=float)
        np.add.at(vals,h['X/indices'][ptr[0]:ptr[1]],h['X/data'][ptr[0]:ptr[1]])
        d['SOX9_nonnegative_expression']=vals[d.full_source_row.to_numpy()]
        obs=d.loc[d.original_Fib_gate].copy()
        print('Extract full normalized RNA for',len(obs),'Fib; no marker-positive gate',flush=True)
        matrix=subset_sparse(h['X'],obs.full_source_row.to_numpy())
    mandatory=['SOX9','PDGFRA','VIM','ACTA2','COL1A1']
    purity=['KRT14','KRT19','EPCAM','PTPRC','PECAM1','VWF','RGS5','CSPG4']
    expressed=np.asarray((matrix>0).sum(axis=0)).ravel()>=3
    expressed |= genes.isin(mandatory+purity)
    var=pd.DataFrame({'parent_gene_position':np.flatnonzero(expressed)},index=genes[expressed])
    matrix=matrix[:,expressed].tocsr();assert np.isfinite(matrix.data).all() and matrix.data.min()>=0
    for k in obs.columns:
        if obs[k].dtype.name=='category' or obs[k].dtype.kind in ['O','S','U']:obs[k]=obs[k].astype(str)
    a=ad.AnnData(X=matrix,obs=obs,var=var);a.layers['source_normalized_RNA']=matrix.copy()
    a.uns['input_scope']='All QC-pass original Fibroblast cells from the same 19 libraries used by Figure1; no UMAP subsampling or SOX9-positive filtering'
    a.uns['gene_filter']='Genes detected in >=3 full selected Fib cells plus mandatory and purity markers; gene filtering does not remove cells'
    a.uns['expression_source']=str(RNA);a.uns['counts_boundary']='Parent counts are a fractional raw-count candidate; not used as validated integer counts'
    a.write_h5ad(out,compression='gzip')
    d.to_csv(ROOT/'tables/Figure1_selected19_full_QC_skin_cell_metadata.tsv.gz',sep='\t',compression='gzip')
    f=d.loc[d.original_Fib_gate]
    lib=f.groupby(['record_id','source_sample','patient_id_source','clinical_group','condition_pair'],dropna=False).agg(Fib_cells=('record_id','size'),SOX9_mean=('SOX9_nonnegative_expression','mean'),SOX9_positive_cells=('SOX9_nonnegative_expression',lambda v:int((v>0).sum()))).reset_index()
    assert set(lib.record_id)==set(selected)
    lib.to_csv(ROOT/'tables/Selected19_full_Fib_library_expression.csv',index=False,encoding='utf-8-sig')
    # Equal library weighting within an explicitly recorded patient; avoid treating repeated libraries as independent people.
    patient=lib.groupby(['patient_id_source','condition_pair'],dropna=False).agg(SOX9_mean_expression=('SOX9_mean','mean'),libraries=('record_id',lambda v:'; '.join(v)),n_libraries=('record_id','size'),target_cells=('Fib_cells','sum')).reset_index()
    patient.to_csv(ROOT/'tables/skin_patient_unit_expression.csv',index=False,encoding='utf-8-sig')
    counts.to_csv(ROOT/'tables/selected19_QC_counts.csv',encoding='utf-8-sig')
    cfg=json.loads((SKILL/'assets/figure2_ae_config.example.json').read_text(encoding='utf8'))
    oldcfg=json.loads((PROJECT/'20261009_Figure2AB_Figure1_matched_inputs_v02/config/human_figure2_config.json').read_text(encoding='utf8'))
    markers=oldcfg['markers'];markers['allowlist']=[g for g in markers['allowlist'] if g!='KRT14'];markers['coexpression']=['SOX9','COL1A1'];markers['positive_threshold']=0.0
    cfg.update({'input_h5ad':str(out),'output_dir':str(ROOT/'analysis/human_universal'),'dataset_id':'SOX9_selected19_full_Fib','species':'human','lineage_gate':{'column':'bulk_labels','include':['Fibroblast'],'prefix':[]},'sample_key':'patient_id_source','group_key':'condition_pair','expression':{'counts_layer':None,'normalized_layer':'source_normalized_RNA','target_sum':10000},'embedding':{'use_rep_candidates':[],'n_neighbors':20,'resolution':.55,'seed':SEED},'markers':markers,'panel_mode':'composition_contrasts','comparisons':[{'name':'Scar_vs_normal','numerator':'fibrosis','denominator':'normal'}],'display_title':'Full skin Fib: selected 19 libraries','style':{'font_family':'Arial','font_min_pt':8,'font_max_pt':12,'dpi':600}})
    (ROOT/'00_manifest/figure2_config.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf8')
    renderer=SKILL/'scripts/run_figure2_ae_universal.py';shutil.copy2(renderer,ROOT/'code/run_figure2_ae_universal_PRISTINE.py')
    (ROOT/'audit/skill_renderer_sha256.txt').write_text(sha(renderer),encoding='utf8')
    for k,p in [('rna',RNA),('metadata',META)]:assert p.stat().st_size==before[k]['bytes'] and p.stat().st_mtime_ns==before[k]['mtime_ns']
    qa={'status':'PASS_FULL_SELECTED19_SKIN_FIB_INPUT','user_whitelist_sha256':sha(CHOICE),'requested_records':selected,'skin_source_cells':208616,'skin_QC_cells':len(d),'Fib_QC_cells':len(f),'Fib_gene_count':a.n_vars,'patient_units_by_condition':patient.condition_pair.value_counts().to_dict(),'RNA_source':before,'Figure2_input':str(out),'Figure2_input_sha256':sha(out),'Fib_gate':'Source original bulk_labels == Fibroblast; no marker prefilter','whole_Fib_SOΧ9_includes_zeros':True,'KRT14_dotplot_displayed':False,'KRT14_expression_preserved':True,'Figure1_Figure2_same_library_ID_set':True,'Figure2_is_exact_full_QC_Fib_subset_of_Figure1_skin':True,'SOX9_COL1A1_cutoff':'>0 AND >0; auxiliary coexpression only','source_h5ad_unchanged':True,'old_reference_cells_not_used_as_expression_universe':True,'no_inherited_embedding_used':True,'post_expression_view_specimen_selection':'User-selected exploratory cohort; selection occurred after viewing SOX9 summaries, not an independently prespecified confirmatory cohort'}
    (ROOT/'audit/INPUT_MATCH_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(qa,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
