from pathlib import Path
import ast, hashlib, json, sys
import h5py, anndata as ad, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_full_selected import ROOT, PROJECT
sys.stdout.reconfigure(encoding='utf8')
ORDER=['skin','lung','kidney','heart','liver']
PALETTE={'Epithelial/parenchymal':'#0072B2','Fibroblast/stromal':'#D55E00','Endothelial':'#009E73','Neural/glial':'#7A7A7A','Myeloid':'#CC79A7','Lymphoid':'#56B4E9','Other':'#C7C7C7'}
PREV=PROJECT/'20261009_Figure1BC_selected7_controls_v02'

def bh(p):
    p=np.asarray(p,float);n=len(p);o=np.argsort(p);q=np.empty(n)
    q[o]=np.minimum(1,np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]);return q

def style():
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})

def main():
    path=ROOT/'02_shared_geometry/Figure1_five_organ_selected19_refit.h5ad'
    with h5py.File(path,'r') as h:
        d=ad.io.read_elem(h['obs']);xy=h['obsm/X_umap'][:]
    assert len(d)==465450 and len(d.loc[d.organ.eq('skin')])==208532
    expected=pd.read_csv(ROOT/'01_input_audit/user_candidate_skin_record_ids_LOCKED.csv').record_id
    assert set(d.loc[d.organ.eq('skin'),'record_id'])==set(expected)
    reference=PROJECT/'20261009_Figure1BC_horizontal_TableS1_v03/code/compose_and_metadata.py'
    text=reference.read_text(encoding='utf8');tree=ast.parse(text)
    funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['reference','raster','show']];assert len(funcs)==3
    limits=np.array([xy.min(axis=0),xy.max(axis=0)])
    ns={'np':np,'pd':pd,'d':d,'xy':xy,'limits':limits,'PALETTE':PALETTE}
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(reference),'exec'),ns)
    d['display_label']=d.apply(ns['reference'],axis=1)
    # Contract appendage lineages for Figure1 display only; preserve original bulk/sub labels.
    appendage=d.organ.eq('skin')&d.bulk_annotation.str.startswith('SAC_')
    d.loc[appendage,'display_label']='Epithelial/parenchymal'
    d[['record_id','organ','cell_annotation','bulk_annotation','display_label']].to_csv(ROOT/'audit/Figure1_display_annotation_ledger.tsv.gz',sep='\t',compression='gzip')
    raster,show=ns['raster'],ns['show']
    previous=pd.read_csv(PREV/'tables/SOX9_target_expression_statistics_exploratory.csv')
    tabs={organ:pd.read_csv(PREV/'tables'/f'{organ}_target_sample_unit_expression.csv') for organ in ORDER if organ!='skin'}
    patient=pd.read_csv(ROOT/'tables/skin_patient_unit_expression.csv',float_precision='round_trip')
    tabs['skin']=patient
    v0=patient.loc[patient.condition_pair.eq('normal'),'SOX9_mean_expression'].to_numpy()
    v1=patient.loc[patient.condition_pair.eq('fibrosis'),'SOX9_mean_expression'].to_numpy()
    assert len(v0)==5 and len(v1)==10 and (v0==0).sum()==2
    p=float(stats.ttest_ind(v0,v1,equal_var=False).pvalue)
    frame=previous.copy();idx=frame.index[frame.organ.eq('skin')][0]
    frame.loc[idx,'P_exploratory']=p;frame.loc[idx,'reference_group']='normal';frame.loc[idx,'case_group']='fibrosis'
    frame.loc[idx,'target']='All Fib'
    frame.loc[idx,'n_reference_sample_units']=len(v0);frame.loc[idx,'n_case_sample_units']=len(v1)
    frame.loc[idx,'mean_reference']=float(v0.mean());frame.loc[idx,'mean_case']=float(v1.mean())
    frame.loc[idx,'delta_case_minus_reference']=float(v1.mean()-v0.mean())
    frame.loc[idx,'test']='Two-sided Welch test on recorded patient-unit means; equal library weighting within patient; user-selected exploratory cohort'
    frame['FDR_exploratory_BH_5']=bh(frame.P_exploratory.to_numpy())
    frame['analysis_units']=frame.organ.map({'skin':'Recorded patient IDs; equal library mean within patient','lung':'Unchanged legacy public sample units','kidney':'Unchanged legacy public sample units','heart':'Unchanged legacy public sample units','liver':'Unchanged legacy public sample units'})
    frame['skin_target_changed_to_all_Fib']=frame.organ.eq('skin')
    frame.to_csv(ROOT/'tables/Figure1C_all_organ_expression_P_FDR.csv',index=False,encoding='utf-8-sig')
    # Verify non-skin inputs and raw P values are untouched. FDR is recomputed because the five-test family changed.
    for organ in ORDER[1:]:
        assert float(frame.loc[frame.organ.eq(organ),'P_exploratory'].iloc[0])==float(previous.loc[previous.organ.eq(organ),'P_exploratory'].iloc[0])
        tabs[organ].to_csv(ROOT/'tables'/f'{organ}_unchanged_public_statistical_units.csv',index=False,encoding='utf-8-sig')
    skin_stats={'target':'All Fib; original Fibroblast annotation; includes zeros','reference_libraries':5,'scar_libraries':14,'reference_patient_IDs_as_recorded':5,'scar_patient_IDs_as_recorded':10,'reference_mean_of_patient_means':float(v0.mean()),'scar_mean_of_patient_means':float(v1.mean()),'mean_difference':float(v1.mean()-v0.mean()),'fold_change_ratio':float(v1.mean()/v0.mean()),'two_sided_Welch_P':p,'BH_FDR_five_organ_family':float(frame.loc[idx,'FDR_exploratory_BH_5']),'two_sided_Mann_Whitney_P_sensitivity':float(stats.mannwhitneyu(v0,v1,alternative='two-sided').pvalue),'post_SOΧ9_review_selected_cohort':'Exploratory; not an independently prespecified confirmatory sample selection','normal_zero_expression_samples_retained':['F091','F092'],'YN3_full_original_Fib_mean':float(patient.loc[patient.libraries.eq('F093'),'SOX9_mean_expression'].iloc[0])}
    skin_stats['Mann_Whitney_sensitivity_BH_FDR_single_test_family']=skin_stats['two_sided_Mann_Whitney_P_sensitivity']
    (ROOT/'audit/Figure1C_SKIN_STATISTICS.json').write_text(json.dumps(skin_stats,ensure_ascii=False,indent=2),encoding='utf8')
    style();W,H=210,75;fig=plt.figure(figsize=(W/25.4,H/25.4))
    fig.text(1/W,73/H,'b',fontsize=12,fontweight='bold',va='top');fig.text(84/W,73/H,'c',fontsize=12,fontweight='bold',va='top')
    ax=fig.add_axes([2/W,9/H,57/W,57/H]);show(ax,raster(1700))
    leg=fig.add_axes([60/W,27/H,24/W,38/H]);leg.set_axis_off();leg.set_xlim(0,1);leg.set_ylim(0,1)
    for i,((label,col),name) in enumerate(zip(PALETTE.items(),['Epi/parenchymal','Fib/stromal','Endothelial','Neural/glial','Myeloid','Lymphoid','Other'])):
        y=.96-i*.15;leg.scatter([.02],[y],s=12,color=col,clip_on=False);leg.text(.105,y,name,fontsize=8,va='center')
    names={'skin':'All Fib','lung':'AT epi','kidney':'Collecting','heart':'Fib','liver':'Epi/Par'}
    limits_y={}
    for i,organ in enumerate(ORDER):
        center=99+i*24;left=center-9
        fig.text(center/W,68/H,organ.title(),fontsize=9,weight='bold',ha='center')
        ax=fig.add_axes([left/W,42/H,18/W,18/H]);show(ax,raster(700,d.organ.eq(organ).to_numpy()))
        fig.text(center/W,35/H,names[organ],fontsize=8,weight='bold',ha='center')
        ax=fig.add_axes([left/W,10/H,18/W,20/H]);row=frame.loc[frame.organ.eq(organ)].iloc[0];tab=tabs[organ]
        vals=[tab.loc[tab.condition_pair.eq(group),'SOX9_mean_expression'].to_numpy() for group in [row.reference_group,row.case_group]]
        assert all(len(v)>0 and np.isfinite(v).all() for v in vals)
        ax.bar([0,1],[v.mean() for v in vals],width=.58,color=['#CBD6E2','#F2E8CA'],edgecolor='black',linewidth=.6)
        rng=np.random.default_rng(20261009+i)
        for j,v in enumerate(vals):ax.scatter(j+rng.uniform(-.1,.1,len(v)),v,s=5,color='#333333',alpha=.7,linewidths=0,zorder=4)
        max_value=max(max(v) for v in vals);raw=max_value*1.45
        tick=10**np.floor(np.log10(max(raw,.01)))/2
        top=np.ceil(raw/tick)*tick;top=max(top,.1);limits_y[organ]=float(top)
        ax.set_ylim(0,top);ax.set_yticks([0,top/2,top],[f'{v:.2g}' for v in [0,top/2,top]]);ax.set_xticks([0,1],['Ref','Fib']);ax.tick_params(length=1.7,pad=1)
        if i==0:ax.set_ylabel(r'$SOX9$ expression',labelpad=2)
        ax.text(.5,.965,f'P={row.P_exploratory:.3g}\nFDR={row.FDR_exploratory_BH_5:.3g}',ha='center',va='top',fontsize=8,transform=ax.transAxes)
    fig.text(4/W,2/H,'Exploratory P/FDR; skin patient units. Organ-specific y-scales; cross-organ UMAP is descriptive.',fontsize=8)
    fig.canvas.draw()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_text() and t.get_visible()]
    sizes=[t.get_fontsize() for t in texts];assert min(sizes)>=8 and max(sizes)<=12
    boxes=[(t.get_text(),t.get_window_extent(fig.canvas.get_renderer())) for t in texts]
    clipped=[s for s,b in boxes if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5]
    stem=ROOT/'figures/Figure1bc_horizontal_selected19_fullcells'
    for ext in ['png','pdf','svg']:fig.savefig(stem.with_suffix('.'+ext),dpi=600,facecolor='white')
    plt.close(fig)
    from PIL import Image
    with Image.open(stem.with_suffix('.png')) as image:image.save(stem.with_suffix('.tiff'),dpi=(600,600),compression='tiff_lzw')
    qa={'status':'RENDERED_REQUIRES_VISUAL_QA','total_UMAP_cells':len(d),'skin_cells':208532,'skin_library_IDs':sorted(expected),'all_cells_used_in_UMAP':True,'new_geometry_not_old':True,'layout_mm':[210,75],'font_pt':[min(sizes),max(sizes)],'clipped_text':clipped,'y_scales':limits_y,'skin_endpoint':'All Fib mean normalized SOX9, zeros included','other_organ_input_cells_expression_units_and_raw_P_unchanged':True,'BH_FDR_recomputed_across_current_five_tests':True,'original_annotation_contraction':{'SAC_family_display':'Epi/parenchymal only; original source labels unchanged'},'source_renderer':str(reference),'source_renderer_sha256':hashlib.sha256(text.encode()).hexdigest(),'annotation_raster_and_show_functions_exact':True}
    (ROOT/'audit/Figure1_RENDER_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(skin_stats,ensure_ascii=False,indent=2),flush=True);print(json.dumps(qa,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
