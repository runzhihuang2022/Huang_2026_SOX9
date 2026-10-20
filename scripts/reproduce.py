#!/usr/bin/env python3
"""Portable reproduction from aggregate inputs. No network, SSH, or Office calls.

Full embedding/clustering is opt-in and requires separately authorized matrices.
The default command never reads a workstation-specific drive or patient key.
"""
from pathlib import Path
import argparse, ast, hashlib, importlib.util, json, os, platform, sys, time
os.environ.setdefault('OMP_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('MKL_NUM_THREADS','4')
os.environ.setdefault('NUMBA_NUM_THREADS','4')
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from PIL import Image

REPO=Path(__file__).resolve().parents[1]
DATA=REPO/'data/processed'
CHECKS=[]
ARTIFACTS=[]
TABLES_ONLY=False

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def read(path):
    return pd.read_csv(path,sep='\t' if Path(path).suffix=='.tsv' else ',',float_precision='round_trip')

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    if isinstance(data,pd.DataFrame):data.to_csv(path,index=False,encoding='utf-8-sig',sep='\t' if path.suffix=='.tsv' else ',')
    else:path.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')

def bh(p):
    p=np.asarray(p,float);q=np.full(p.shape,np.nan);ok=np.isfinite(p);v=p[ok];n=len(v)
    if n:
        o=np.argsort(v);z=np.empty(n);z[o]=np.minimum(1,np.minimum.accumulate((v[o]*n/np.arange(1,n+1))[::-1])[::-1]);q[ok]=z
    return q

def check(name,actual,expected,rtol=1e-10,atol=1e-12):
    good=bool(np.allclose(actual,expected,rtol=rtol,atol=atol,equal_nan=True))
    CHECKS.append({'check':name,'passed':good,'rtol':rtol,'atol':atol})
    if not good:raise AssertionError(name)

def style():
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
        'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':.6,
        'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom',
        'mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})

def save(fig,out,name):
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    sizes=[t.get_fontsize() for t in texts]
    clipped=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if b.x0<-1 or b.y0<-1 or b.x1>fig.bbox.width+1 or b.y1>fig.bbox.height+1:clipped.append(t.get_text())
    # Historical Figure5 may clip labels; preserve and report rather than silently redesign it.
    for ext in ['pdf','svg','png']:fig.savefig(out/'figures'/f'{name}.{ext}',dpi=600,facecolor='white')
    ARTIFACTS.append({'asset':name,'canvas_mm':(fig.get_size_inches()*25.4).tolist(),
        'font_pt':[min(sizes),max(sizes)] if sizes else [],'clipped_text':clipped,'dpi':600})
    plt.close(fig)

def legacy_ast(relative):
    p=REPO/relative;return p,ast.parse(p.read_text(encoding='utf-8'))

def figure1(out):
    order=['skin','lung','kidney','heart','liver'];expected=read(DATA/'figure1/expected_statistics.csv')
    frame=expected.copy();tabs={organ:read(DATA/f'figure1/{organ}_statistical_units.csv') for organ in order}
    for i,r in frame.iterrows():
        d=tabs[r.organ];x=d.loc[d.condition_pair.eq(r.reference_group),'SOX9_mean_expression'].to_numpy()
        y=d.loc[d.condition_pair.eq(r.case_group),'SOX9_mean_expression'].to_numpy()
        check(r.organ+'_n',[len(x),len(y)],[r.n_reference_sample_units,r.n_case_sample_units],0,0)
        check(r.organ+'_Welch_P',stats.ttest_ind(x,y,equal_var=False).pvalue,r.P_exploratory)
        check(r.organ+'_group_means',[x.mean(),y.mean()],[r.mean_reference,r.mean_case])
        frame.loc[i,'P_exploratory']=float(stats.ttest_ind(x,y,equal_var=False).pvalue)
    frame['FDR_exploratory_BH_5']=bh(frame.P_exploratory)
    check('five_organ_BH',frame.FDR_exploratory_BH_5,expected.FDR_exploratory_BH_5)
    check('retained_normal_zero_expression_units',sum(tabs['skin'].SOX9_mean_expression.eq(0)),2,0,0)
    write(out/'tables/Figure1C_recomputed_statistics.csv',frame)
    if TABLES_ONLY:return frame
    # Execute the approved plotting block unchanged, with frozen density rasters as display input.
    p,tree=legacy_ast('scripts/legacy/atlas/render_Figure1BC.py')
    main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and getattr(n.value.func,'id',None)=='style')
    stop=next(i+1 for i,n in enumerate(main.body[start:],start) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and n.value.func.attr=='close')
    def raster(size,focus=None):
        organ='all' if focus is None else order[int(np.flatnonzero(focus)[0])]
        with Image.open(REPO/f'data/frozen_maps/figure1_{organ}.png') as im:return np.asarray(im)
    def show(ax,img):ax.imshow(img,interpolation='nearest',rasterized=True);ax.set_axis_off();ax.set_aspect('equal')
    palette={'Epithelial/parenchymal':'#0072B2','Fibroblast/stromal':'#D55E00','Endothelial':'#009E73','Neural/glial':'#7A7A7A','Myeloid':'#CC79A7','Lymphoid':'#56B4E9','Other':'#C7C7C7'}
    ns={'np':np,'pd':pd,'plt':plt,'matplotlib':matplotlib,'ROOT':out,'ORDER':order,
        'PALETTE':palette,'frame':frame,'tabs':tabs,'style':style,'raster':raster,'show':show,
        'd':pd.DataFrame({'organ':order})}
    exec(compile(ast.Module(body=main.body[start:stop],type_ignores=[]),str(p),'exec'),ns)
    ARTIFACTS.append({'asset':'Figure1bc_horizontal_selected19_fullcells','mode':'Original plotting block; density-map display frozen, Welch/BH recomputed'})
    return frame

def figure2(out):
    marker=read(DATA/'figure2/Figure2_marker_allowlist_summary.tsv')
    ledger=read(DATA/'figure2/Figure2_selected_marker_ledger.tsv')
    info=json.loads((REPO/'config/figure2_display.json').read_text())
    clusters=info['clusters'];mandatory=['SOX9','PDGFRA','VIM','ACTA2']
    displayed=list(dict.fromkeys(mandatory+ledger.gene.tolist()))
    assert 'KRT14' not in displayed
    purity=read(DATA/'figure2/Figure2_all_cluster_marker_and_purity_audit.csv')
    sox=purity[purity.gene.eq('SOX9')].sort_values('mean_expression',ascending=False)
    assert sox.iloc[0].cluster=='C12'
    c12=purity[purity.cluster.eq('C12') & purity.gene.isin(['SOX9','ACTA2','VIM'])].set_index('gene')
    check('C12_gene_detection_percent',c12.loc[['SOX9','ACTA2','VIM'],'percent_expressed'],np.array([61,180,178])/185*100)
    check('Figure2_total_Fib',purity[purity.gene.eq('SOX9')].n_cells.sum(),20270,0,0)
    write(out/'tables/Figure2_C12_marker_percentages.csv',c12.reset_index())
    # Recompute the registered sample-level composition contrasts (all 17 tests).
    summary=read(DATA/'figure2/Figure2_sample_cluster_summary.tsv')
    p,tree=legacy_ast('scripts/renderers/run_figure2_ae_universal.py')
    funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['bh','contrast_table','highest_cluster_tests']]
    ns={'pd':pd,'np':np,'stats':stats};exec(compile(ast.Module(body=funcs,type_ignores=[]),str(p),'exec'),ns)
    # Preserve the pristine routine's complete sample x cluster grid, including zero cells.
    tests=ns['contrast_table'](summary,{'name':'Scar_vs_normal','numerator':'fibrosis','denominator':'normal'})
    expected=read(DATA/'figure2/Figure2D_Scar_vs_normal_P_FDR.tsv').set_index('cluster').loc[clusters]
    for c in ['Mann_Whitney_P','BH_FDR','mean_fraction_difference']:check('Figure2D_'+c,tests[c],expected[c])
    write(out/'tables/Figure2D_composition_recomputed.csv',tests)
    high=ns['highest_cluster_tests'](summary,'C12')
    expected_high=read(DATA/'figure2/Figure2_highest_target_cluster_tests_P_FDR.tsv')
    for c in ['mean_paired_difference','one_sided_paired_Wilcoxon_P','BH_FDR']:check('Figure2_highest_'+c,high[c],expected_high[c])
    write(out/'tables/Figure2_highest_cluster_recomputed.csv',high)
    if TABLES_ONLY:return
    style();W,H=156.66,60.0
    # Use the stored actual panel dimensions, not a new arbitrary canvas.
    import fitz
    doc=fitz.open(REPO/'outputs/reference/Figure2AB_human_Fib_slot.pdf');W,H=doc[0].rect.width/72*25.4,doc[0].rect.height/72*25.4;doc.close()
    fig=plt.figure(figsize=(W/25.4,H/25.4))
    fig.text(1/W,(H-1)/H,'a',fontsize=12,weight='bold',va='top');fig.text(61/W,(H-1)/H,'b',fontsize=12,weight='bold',va='top')
    fig.text(8/W,(H-1)/H,'Human Fib clusters',fontsize=9,weight='bold',va='top');fig.text(68/W,(H-1)/H,'Marker expression',fontsize=9,weight='bold',va='top')
    ax=fig.add_axes([3/W,23/H,48/W,33/H]);ax.imshow(Image.open(REPO/'data/frozen_maps/figure2_clusters.png'));ax.axis('off')
    handles=[Line2D([],[],marker='o',ls='',ms=2,color=info['palette'][c],label=c) for c in clusters]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.07,-.02),ncol=5,frameon=False,fontsize=8,
              markerscale=3,handlelength=.8,handletextpad=.25,columnspacing=.55,borderaxespad=0)
    fig.text(4/W,1.2/H,info['highest_cluster']+r'*: highest mean $SOX9$',fontsize=8,va='bottom')
    dot_axis=fig.add_axes([77/W,11/H,67/W,40/H])
    p,tree=legacy_ast('scripts/renderers/run_figure2_ae_universal.py')
    main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='dot' for t in n.targets))
    stop=next(i for i,n in enumerate(main.body[start:],start) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='panel_contract' for t in n.targets))
    ns={'fig':fig,'dot_axis':dot_axis,'dot_letter':'b','marker':marker,'displayed_genes':displayed,'cluster_order':clusters,'mandatory':mandatory,'np':np}
    exec(compile(ast.Module(body=main.body[start:stop],type_ignores=[]),str(p),'exec'),ns)
    dot_axis.set_title('',loc='left');dot_axis.set_xlabel('')
    for t in dot_axis.get_xticklabels():t.set_rotation(90);t.set_horizontalalignment('center')
    cax=fig.add_axes([78/W,55/H,27/W,1.8/H]);cb=fig.colorbar(ScalarMappable(norm=Normalize(-2,2),cmap='RdBu_r'),cax=cax,orientation='horizontal')
    cb.set_ticks([-2,0,2]);cb.ax.tick_params(labelsize=8,length=1,pad=.5)
    fig.text(109/W,56/H,'Mean expression (z)',fontsize=8,va='center');fig.text(77/W,1/H,'Dot size: percent expressing',fontsize=8,va='bottom')
    save(fig,out,'Figure2AB_from_aggregate_markers_frozen_map')

def spatial(out):
    h=read(DATA/'spatial/Human_marker_by_source_part_region.csv')
    merged=h.groupby(['section','region'],sort=False)[['n_source_Fib','dual_positive']].sum().reset_index()
    merged['percent']=merged.dual_positive/merged.n_source_Fib*100;merged['P_value']='NE';merged['BH_FDR']='NE'
    expected=read(DATA/'spatial/Human_marker_by_physical_section_region.csv')
    a=merged.set_index(['section','region']).sort_index();b=expected.set_index(['section','region']).sort_index()
    check('human_region_counts',a[['n_source_Fib','dual_positive','percent']],b[['n_source_Fib','dual_positive','percent']])
    r=read(DATA/'spatial/Figure3G_main_sections_counts.csv')
    check('rat_dual_counts',r.double_positive_Fib,[52,304,160,24,30,28],0,0)
    r['dual_percent']=100*r.double_positive_Fib/r.source_Fib_observations
    pooled=r.groupby('stage')[['source_Fib_observations','double_positive_Fib']].sum().reset_index()
    pooled['dual_percent']=100*pooled.double_positive_Fib/pooled.source_Fib_observations
    pooled['P_value']='NE';pooled['BH_FDR']='NE'
    write(out/'tables/Human_regions_recomputed.csv',merged);write(out/'tables/Rat_sections_recomputed.csv',r)
    write(out/'tables/Rat_stage_pooled_descriptive.csv',pooled)
    if TABLES_ONLY:return
    # Main spatial maps are byte-identical approved reference panels, not recalculated from PNGs.
    # Their external-data regeneration is documented in docs/FULL_ANALYSIS.md.
    style();fig,axes=plt.subplots(1,2,figsize=(210/25.4,70/25.4));fig.subplots_adjust(left=.08,right=.97,bottom=.27,top=.85,wspace=.4)
    for i,section in enumerate(['2 mph','6 mph']):
        vals=merged[merged.section.eq(section)].set_index('region').loc[['Scarless side','Scar side'],'percent']
        axes[0].plot([0,1],vals,marker='o' if i==0 else 's',label=section,color=['#8F7A16','#B29A37'][i])
    axes[0].set_xticks([0,1],['Scarless side','Scar side']);axes[0].set_ylabel('Dual-positive Fib (%)');axes[0].legend(frameon=False)
    for i,stage in enumerate(['Adult','E16.5']):
        d=r[r.stage.eq(stage)];axes[1].scatter(np.full(len(d),i),d.dual_percent,s=16,c='#E6B800',edgecolors='#555555',lw=.4)
        axes[1].plot([i-.15,i+.15],[d.dual_percent.mean()]*2,color='#555555')
    axes[1].set_xticks([0,1],['Adult','E16.5']);axes[1].set_ylabel('Dual-positive Fib (%)')
    for ax in axes:ax.set_title('Descriptive; P=NE, FDR=NE',fontsize=8);ax.spines[['top','right']].set_visible(False)
    fig.text(.04,.04,'Spatial observations are not independent donors/animals. Section times differ; no group hypothesis test.',fontsize=8)
    save(fig,out,'Spatial_counts_validation_not_main_panel')

def culture(out):
    pb=read(DATA/'culture/Figure5A_library_state_pseudobulk_log2CPM.tsv')
    genes=['SCD','TFEB','CTSB','MMP14','MRC2'];rows=[]
    check('culture_libraries',pb.library.nunique(),8,0,0)
    for state,q in pb.groupby('state',sort=False):
        for g in genes:
            rho,p=stats.spearmanr(q.SOX9,q[g]);rows.append({'state':state,'gene':g,'spearman_rho':rho,'P_value':p})
    cors=pd.DataFrame(rows);cors['BH_FDR_within_state']=cors.groupby('state').P_value.transform(lambda v:bh(v))
    expected=read(DATA/'culture/Figure5A_SOX9_target_library_correlations.tsv').set_index(['state','gene']).sort_index()
    actual=cors.set_index(['state','gene']).sort_index()
    for c in ['spearman_rho','P_value','BH_FDR_within_state']:check('culture_correlations_'+c,actual[c],expected[c])
    rows=[]
    for g in ['SOX9']+genes:
        w=pb.pivot(index='library',columns='state',values=g).dropna();delta=w.Fib_ACTA2-w['Other Fib']
        rows.append({'gene':g,'paired_mean_difference':delta.mean(),'paired_median_difference':delta.median(),
            'Fib_ACTA2_mean_log2CPM':w.Fib_ACTA2.mean(),'Other_Fib_mean_log2CPM':w['Other Fib'].mean(),
            'paired_exact_Wilcoxon_P':stats.wilcoxon(w.Fib_ACTA2,w['Other Fib'],alternative='two-sided',method='exact').pvalue})
    diffs=pd.DataFrame(rows);diffs['BH_FDR']=bh(diffs.paired_exact_Wilcoxon_P)
    expected=read(DATA/'culture/Figure5A_Fib_ACTA2_vs_other_Fib_paired_differences.tsv').set_index('gene').loc[diffs.gene]
    for c in ['paired_mean_difference','paired_exact_Wilcoxon_P','BH_FDR']:check('culture_paired_'+c,diffs[c],expected[c])
    write(out/'tables/Figure5A_correlations_recomputed.csv',cors);write(out/'tables/Figure5A_paired_differences_recomputed.csv',diffs)
    if TABLES_ONLY:return
    # Original Figure5 rendering block, no inferred replacement or figure digitization.
    p,tree=legacy_ast('scripts/legacy/culture/run_server_figure4h_5a.py')
    main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
               and isinstance(n.value.func,ast.Attribute) and ast.unparse(n.value.func)=='mpl.rcParams.update')
    stop=next(i+1 for i,n in enumerate(main.body[start:],start) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
              and isinstance(n.value.func,ast.Attribute) and ast.unparse(n.value.func)=='plt.close')
    ns={'np':np,'pd':pd,'mpl':matplotlib,'plt':plt,'pb':pb,'cors':cors,'diffs':diffs,'GENES_5A':genes,'fig_dir':out/'figures'}
    exec(compile(ast.Module(body=main.body[start:stop],type_ignores=[]),str(p),'exec'),ns)
    ARTIFACTS.append({'asset':'Figure5A_SOX9_SCD_lysosomal_targets_correlations_and_differences_A4width','mode':'Original plotting block; eight-library correlations and paired tests recomputed'})
    spec=importlib.util.spec_from_file_location('net',REPO/'scripts/renderers/sox9_network_go.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.configure();m.SQUARE=(105*m.MM,105*m.MM)
    edges=read(DATA/'culture/Figure4H_Fib_ACTA2_SOX9_motif_pruned_targets_with_correlations.tsv')
    go=read(DATA/'culture/Figure4H_Fib_ACTA2_SOX9_regulon_GO_BP.tsv')
    check('SOX9_target_edges',len(edges),37,0,0)
    for label,fig in [('Figure4_SOX9_target_network',m.network_figure(edges,'Fib_ACTA2','SOX9')),
                       ('Figure4_SOX9_GO_cached',m.go_figure(go,'Fib_ACTA2','SOX9')[0])]:save(fig,out,label)
    # Cached enrichment cannot certify missing upstream pySCENIC execution or GO membership.

def full_figure2(args,out):
    if not args.input:raise ValueError('--input authorized H5AD is required')
    import subprocess
    cfg=json.loads((REPO/'config/figure2_full_analysis.json').read_text())
    cfg['input_h5ad']=str(Path(args.input).resolve());cfg['output_dir']=str(out/'full_figure2')
    cfg['sample_key']=args.sample_key
    path=out/'figure2_effective_config.json';write(path,cfg)
    subprocess.run([sys.executable,str(REPO/'scripts/renderers/run_figure2_ae_universal.py'),'--config',str(path)],check=True)
    actual=json.loads((out/'full_figure2/00_manifest/Figure2_A-E_manifest.json').read_text())
    expected=json.loads((REPO/'config/figure2_display.json').read_text())
    check('full_Figure2_cell_gene_shape',actual['lineage_shape'],[20270,35379],0,0)
    check('full_Figure2_cluster_count',actual['n_clusters'],17,0,0)
    if actual['highest_target_cluster']!=expected['highest_cluster']:raise AssertionError('Highest cluster differs')
    CHECKS.append({'check':'full_Figure2_highest_C12','passed':True})
    # Anonymous barcode replacement changes a barcode-table hash even with identical coordinates.
    ARTIFACTS.append({'asset':'full_Figure2_geometry','original_barcode_table_hash_equal':actual['geometry_sha256']==expected['geometry_sha256_original_barcode_table'],
        'interpretation':'Anonymous barcode replacement changes this hash; coordinate-array fingerprint is tested separately'})
    # gzip is inferred by pandas; this table is tab-separated, not CSV.
    geo=pd.read_csv(out/'full_figure2/02_shared_geometry/SOX9_selected19_full_Fib_lineage_umap.tsv.gz',sep='\t',float_precision='round_trip')
    numeric_hash=hashlib.sha256(geo[['UMAP1','UMAP2']].to_numpy(dtype='<f8').tobytes()).hexdigest()
    CHECKS.append({'check':'full_Figure2_numeric_geometry_hash','passed':numeric_hash==expected['numeric_geometry_sha256_float64_ordered']})
    current=read(out/'full_figure2/05_tables/Figure2_marker_allowlist_summary.tsv').set_index(['cluster','gene']).sort_index()
    reference=read(DATA/'figure2/Figure2_marker_allowlist_summary.tsv').set_index(['cluster','gene']).sort_index()
    for c in ['mean_expression','percent_expressed','n_cells']:check('full_Figure2_'+c,current[c],reference[c],rtol=1e-7,atol=1e-10)

def full_atlas(args,out):
    if not args.input:raise ValueError('--input authorized 465450 x 1679 normalized analysis H5AD required')
    import anndata as ad,scanpy as sc
    a=ad.read_h5ad(args.input);assert a.shape==(465450,1679)
    # Drop prior graph/embedding so a saved embedding is never reused as a fresh analysis.
    a.obsm.clear();a.obsp.clear()
    sc.pp.pca(a,n_comps=50,zero_center=True,svd_solver='arpack',random_state=20261009)
    sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca',random_state=20261009,transformer='pynndescent')
    sc.tl.umap(a,random_state=20261009,min_dist=.35,spread=1.0,maxiter=200)
    sc.tl.leiden(a,resolution=.5,key_added='leiden_global_new',random_state=20261009,flavor='igraph',n_iterations=2,directed=False)
    a.write_h5ad(out/'full_atlas.h5ad',compression='gzip')

def full_culture(args,out):
    if not args.input:raise ValueError('--input authorized cultured-Fib raw-count H5AD required')
    p=REPO/'scripts/legacy/culture/run_server_figure4h_5a.py'
    spec=importlib.util.spec_from_file_location('full_culture',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.OUT=out;m.H5AD=Path(args.input).resolve()
    m.REGULON=REPO/'data/processed/culture/s3_avg20_rep1.regulons.txt'
    m.RENDERER=REPO/'scripts/renderers/sox9_network_go.py'
    # Freeze the original query result; a live service update must not silently alter the paper.
    m.gprofiler_go=lambda genes,background:read(DATA/'culture/Figure4H_Fib_ACTA2_SOX9_regulon_GO_BP.tsv')
    m.main()
    for name,keys in [('Figure5A_library_state_pseudobulk_log2CPM.tsv',['library','state']),
                      ('Figure4H_Fib_ACTA2_SOX9_motif_pruned_targets_with_correlations.tsv',['target_symbol'])]:
        actual=read(out/'05_tables'/name).set_index(keys).sort_index()
        expected=read(DATA/'culture'/name).set_index(keys).sort_index()
        for c in expected.select_dtypes(include='number').columns:check('full_culture_'+name+'_'+c,actual[c],expected[c],rtol=1e-8,atol=1e-10)

def main():
    global TABLES_ONLY
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task',choices=['all','figure1','figure2','spatial','culture','full-figure2','full-atlas','full-culture'])
    parser.add_argument('--out',required=True,help='A NEW output directory; existing directories are refused')
    parser.add_argument('--input',help='Authorized external H5AD for opt-in full mode')
    parser.add_argument('--sample-key',default='public_subject_id',help='Subject key in the authorized full Figure2 input')
    parser.add_argument('--tables-only',action='store_true',help='Validate numbers without plotting; useful on CI without licensed Arial')
    args=parser.parse_args();out=Path(args.out).resolve()
    TABLES_ONLY=args.tables_only
    if not TABLES_ONLY:
        from matplotlib.font_manager import FontProperties,findfont
        findfont(FontProperties(family='Arial'),fallback_to_default=False)
    if out.exists():raise FileExistsError('Refusing to overwrite '+str(out))
    out.mkdir(parents=True);(out/'tables').mkdir();(out/'figures').mkdir()
    started=time.time();errors=[]
    for name,func in [('figure1',figure1),('figure2',figure2),('spatial',spatial),('culture',culture)]:
        if args.task in ['all',name]:
            print('RUN',name,flush=True)
            try:func(out)
            except Exception as e:
                import traceback;traceback.print_exc();errors.append({'task':name,'error':str(e)})
    if args.task=='full-figure2':full_figure2(args,out)
    if args.task=='full-atlas':full_atlas(args,out)
    if args.task=='full-culture':full_culture(args,out)
    checksums=[{'path':str(p.relative_to(out)).replace('\\','/'),'sha256':sha(p),'bytes':p.stat().st_size} for p in out.rglob('*') if p.is_file()]
    passed=not errors and all(c['passed'] for c in CHECKS)
    result={'status':('executed_verified' if CHECKS else 'partial') if passed else 'executed_failed',
        'scope':'Aggregate numerical replay and documented frozen-display rendering; NOT raw-data-to-paper verification',
        'task':args.task,'checks':CHECKS,'errors':errors,'figure_checks':ARTIFACTS,'files':checksums,
        'elapsed_seconds':time.time()-started,'python':platform.python_version(),'platform':platform.platform(),
        'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}
    write(out/'reproduction_report.json',result)
    print(json.dumps({'status':result['status'],'numerical_checks':len(CHECKS),'errors':errors,'elapsed_seconds':result['elapsed_seconds']},ensure_ascii=False))
    if result['status']=='executed_failed':sys.exit(1)

if __name__=='__main__':main()
