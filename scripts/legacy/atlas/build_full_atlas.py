from pathlib import Path
import gc, json, hashlib, sys, time
import numpy as np, pandas as pd, h5py, anndata as ad, scanpy as sc
from scipy import sparse
from prepare_full_selected import ROOT, PROJECT, RNA, subset_sparse, sha
sys.stdout.reconfigure(encoding='utf8')
SEED=20261009
AUDIT=PROJECT/'20261009_Figure1ABC_all_samples_Fib_SOX9_v01/audit/live_source_audit.json'
PREV=PROJECT/'20261009_Figure1BC_selected7_controls_v02'

def symbols(v):
    s=v['gene_symbol_original'].astype(str).str.upper() if 'gene_symbol_original' in v else pd.Series(v.index.astype(str).str.upper(),index=v.index)
    s=s[~s.isin(['','NAN','NONE'])];return s

def main():
    out=ROOT/'02_shared_geometry/Figure1_five_organ_selected19_refit.h5ad'
    if out.exists():raise FileExistsError('Completed atlas exists; do not overwrite')
    audit=json.loads(AUDIT.read_text(encoding='utf8'))
    skin=pd.read_csv(ROOT/'tables/Figure1_selected19_full_QC_skin_cell_metadata.tsv.gz',sep='\t',index_col=0,low_memory=False,float_precision='round_trip')
    old=pd.read_csv(PREV/'audit/selected_cells_metadata_expression.tsv.gz',sep='\t',index_col=0,low_memory=False,float_precision='round_trip')
    public=old.loc[~old.organ.eq('skin')].copy();assert len(public)==256918
    registries={};inputs={};before={}
    for k in [2,3,4,5,6]:
        p=Path(audit[str(k)]['path'].replace('prepared_hvg_cluster_annotated','combined_raw_counts'));inputs[k]=p
        before[str(k)]={'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
        with h5py.File(p,'r') as h:v=ad.io.read_elem(h['var'])
        s=symbols(v);unique=s[~s.duplicated(keep=False)]
        registries[k]=pd.Series(np.arange(len(v)),index=v.index).reindex(unique.index).set_axis(unique.to_numpy())
    with h5py.File(RNA,'r') as h:
        vg=ad.io.read_elem(h['var']);human=pd.Index(vg.index.astype(str))
    common=set(human)
    for registry in registries.values():common &= set(registry.index)
    assert len(common)>=5000,len(common)
    features=[];ledger=[]
    for k in [0,2,3,4,5,6]:
        with h5py.File(audit[str(k)]['path'],'r') as h:v=ad.io.read_elem(h['var'])
        v=v.copy();v['symbol']=v.index.astype(str).str.upper();v=v.loc[v.symbol.isin(common)]
        rank=next((c for c in ['highly_variable_rank','dispersions_norm','variances_norm'] if c in v.columns),None)
        if rank:
            v=v.sort_values(rank,ascending=rank=='highly_variable_rank',na_position='last',kind='stable')
        # A fixed, balanced feature union; no cell sampling and no SOX9-effect based tuning.
        selected=v.symbol.drop_duplicates().head(500).tolist()
        features.extend(selected)
        ledger.extend({'source_index':k,'symbol':g,'within_source_feature_rank':i+1,'source_rank_field':rank or 'existing source HVG order'} for i,g in enumerate(selected))
    markers=['SOX9','PDGFRA','VIM','ACTA2','COL1A1','COL3A1','DCN','LUM','KRT14','KRT19','EPCAM','PTPRC','PECAM1','VWF','RGS5']
    features=list(dict.fromkeys(features+[g for g in markers if g in common]));assert len(features)>1200
    pd.DataFrame(ledger).to_csv(ROOT/'audit/Atlas_feature_selection_ledger.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame({'symbol':features}).to_csv(ROOT/'audit/Atlas_shared_expression_features.csv',index=False,encoding='utf-8-sig')
    matrices=[];observations=[]
    print('Atlas features',len(features),'full shared unique symbol pool',len(common),flush=True)
    with h5py.File(RNA,'r') as h:
        pos=human.get_indexer(features);assert (pos>=0).all()
        order=np.argsort(pos);m=subset_sparse(h['X'],skin.full_source_row.to_numpy(),pos[order])[:,np.argsort(order)].tocsr()
    # SOX9 and every other expression value retain the full-source normalized RNA scale.
    if 'SOX9' in features:assert np.allclose(m[:,features.index('SOX9')].toarray().ravel(),skin.SOX9_nonnegative_expression,rtol=1e-6,atol=1e-6)
    skin.index=skin.index.astype(str)+'|skin_full';matrices.append(m);observations.append(skin)
    print('Atlas skin rows',len(skin),flush=True)
    for k in [2,3,4,5,6]:
        d=public.loc[public.source_index.eq(k)].copy()
        assert len(d)>0
        with h5py.File(inputs[k],'r') as h:
            rawobs=ad.io.read_elem(h['obs']);idx=pd.Index(rawobs.index.astype(str));barcodes=d.index.str.replace(r'-'+str(k)+r'$','',regex=True)
            rows=idx.get_indexer(barcodes);assert (rows>=0).all() and len(np.unique(rows))==len(rows)
            cols=registries[k].reindex(features).to_numpy(dtype=int)
            order=np.argsort(cols);counts=subset_sparse(h['X'],rows,cols[order])[:,np.argsort(order)].tocsr()
            assert counts.data.min()>=0
            denominator=rawobs.iloc[rows]['n_counts'].to_numpy(dtype=float);assert (denominator>0).all()
            x=counts.multiply((10000/denominator).astype(np.float32)[:,None]).tocsr();x.data=np.log1p(x.data)
            if 'SOX9' in features:
                error=np.max(np.abs(x[:,features.index('SOX9')].toarray().ravel()-d.SOX9_nonnegative_expression.to_numpy()))
                assert error<1e-5,(k,error)
            d['raw_parent_row']=rows
        matrices.append(x);observations.append(d)
        del counts;gc.collect()
        print('Atlas public source',k,'rows',len(d),flush=True)
    X=sparse.vstack(matrices,format='csr',dtype=np.float32)
    obs=pd.concat(observations,axis=0,sort=False)
    assert len(obs)==465450 and obs.index.is_unique
    assert set(obs.loc[obs.organ.eq('skin'),'record_id'])==set(skin.record_id)
    # Stringify only metadata to avoid unsupported mixed categoricals; numeric expression remains numeric.
    for c in obs.columns:
        if obs[c].dtype.kind in ['O','S','U'] or obs[c].dtype.name=='category':obs[c]=obs[c].fillna('NR').astype(str)
    a=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=features))
    del matrices,observations,X,m,x;gc.collect()
    a.uns['scope_warning']='Cross-organ human/mouse visualization is descriptive; exact source-symbol uppercase matching is not a validated one-to-one orthology map. Source/statistical means are not directly compared between organs.'
    a.uns['integration_method']='De-novo joint centered PCA on unmodified normalized expression; no batch correction of the organ/species-confounded dataset key.'
    a.uns['feature_rule']='Union of up to 500 pre-existing within-source HVGs per source in shared unique-symbol pool plus registered lineage markers; all chosen QC cells retained'
    a.uns['old_global_X_umap_used']=False;a.uns['all_cells_used_in_embedding']=True
    print('PCA on ALL',a.shape,'joint normalized RNA rows',flush=True)
    sc.pp.pca(a,n_comps=50,zero_center=True,svd_solver='arpack',random_state=SEED)
    print('PCA complete; neighbors on ALL cells',flush=True)
    np.save(ROOT/'02_shared_geometry/Figure1_all_cells_joint_PCA.npy',a.obsm['X_pca'])
    a.obs.to_csv(ROOT/'tables/Figure1_selected19_plus_unchanged_public_cell_metadata.tsv.gz',sep='\t',compression='gzip')
    sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca',random_state=SEED,transformer='pynndescent')
    print('Neighbor graph complete; fresh UMAP on ALL cells',flush=True)
    sc.tl.umap(a,random_state=SEED,min_dist=.35,spread=1.0,maxiter=200)
    sc.tl.leiden(a,resolution=.5,key_added='leiden_global_new',random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
    a.write_h5ad(out,compression='gzip')
    g=pd.DataFrame({'barcode':a.obs_names,'UMAP1':a.obsm['X_umap'][:,0],'UMAP2':a.obsm['X_umap'][:,1],'record_id':a.obs.record_id.to_numpy(),'organ':a.obs.organ.to_numpy()})
    g.to_csv(ROOT/'02_shared_geometry/Figure1_five_organ_full_cell_UMAP.tsv.gz',sep='\t',index=False,compression='gzip')
    ghash=hashlib.sha256(g[['barcode','UMAP1','UMAP2']].to_csv(sep='\t',index=False).encode()).hexdigest()
    for k,p in inputs.items():assert p.stat().st_size==before[str(k)]['bytes'] and p.stat().st_mtime_ns==before[str(k)]['mtime_ns']
    qa={'status':'PASS_FULL_CELL_DE_NOVO_UMAP','shape':list(a.shape),'skin_cells':len(skin),'public_cells_unchanged':len(public),'total_cells':len(a),'record_count':int(a.obs.record_id.nunique()),'skin_library_count':19,'SOX9_gene_expression_not_embedding':True,'old_global_X_umap_used':False,'all_cells_used_no_cell_sampling':True,'source_files_unchanged':True,'source_inputs':before,'geometry_sha256':ghash,'seed':SEED,'n_neighbors':20,'UMAP_epochs':200,'PCA_components':50,'method':'Uncorrected joint normalized RNA PCA; dataset/organ/species are confounded, not corrected away','cross_species_boundary':str(a.uns['scope_warning']),'gene_matching':'Source symbol uppercase; duplicate symbols excluded from embedding features','shared_unique_symbols':len(common)}
    (ROOT/'audit/Figure1_GLOBAL_UMAP_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(qa,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
