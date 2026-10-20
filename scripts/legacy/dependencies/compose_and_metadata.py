from pathlib import Path
import json, hashlib, html, sys
import numpy as np, pandas as pd, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf8')
OUT=Path(__file__).resolve().parents[1]
PREV=Path('data/external/INPUT_8921332a6e/20261009_Figure1BC_selected7_controls_v02')
ORIGINAL=Path('data/external/INPUT_bac603486c/20261009_Figure1ABC_all_samples_Fib_SOX9_v01')
for folder in ['audit','tables','figures','report','outputs']:(OUT/folder).mkdir(parents=True,exist_ok=True)
ORDER=['skin','lung','kidney','heart','liver']
PALETTE={'Epithelial/parenchymal':'#0072B2','Fibroblast/stromal':'#D55E00','Endothelial':'#009E73','Neural/glial':'#7A7A7A','Myeloid':'#CC79A7','Lymphoid':'#56B4E9','Other':'#C7C7C7'}
TARGETS={'skin':['Fib_Myo'],'lung':['Alveolar epithelial cells'],'kidney':['Collecting duct epithelial cells'],'heart':['Fibroblasts','Activated fibroblasts/myofibroblasts']}
NAMES={'skin':'Fib_Myo','lung':'AT epi','kidney':'Collecting','heart':'Fib','liver':'Epi/Par'}
AUD=json.loads((ORIGINAL/'audit/live_source_audit.json').read_text(encoding='utf8'))
for entry in AUD.values():
    s=Path(entry['path']).stat();assert s.st_size==entry['bytes'] and s.st_mtime_ns==entry['mtime_ns']
d=pd.read_csv(PREV/'audit/selected_cells_metadata_expression.tsv.gz',sep='\t',index_col=0,low_memory=False)
m=pd.read_csv(PREV/'tables/included_source_sample_metadata.csv')
assert len(d)==312292 and len(m)==144 and set(d.record_id)==set(m.record_id)
assert set(d.loc[(d.organ=='skin')&(d.condition_pair=='normal'),'record_id'])=={'F086','F087','F089','F090','F091','F092','F093'}
with h5py.File(AUD['atlas']['path'],'r') as h:
    index_group=h['obs'][h['obs'].attrs['_index']]
    idx=index_group.asstr()[:] if isinstance(index_group,h5py.Dataset) else index_group['values'].asstr()[:]
    if isinstance(index_group,h5py.Group):assert not index_group['mask'][:].any()
    fullxy=h['obsm/X_umap'][:]
    loc=pd.Index(idx).get_indexer(d.index);assert np.all(loc>=0);xy=fullxy[loc]
limits=np.array([fullxy.min(axis=0),fullxy.max(axis=0)])

def reference(row):
    text=(str(row.bulk_annotation)+' '+str(row.cell_annotation)).lower()
    rules=[('Fibroblast/stromal',['fib','stromal','myofib','stellate','pericyte','smooth']),('Endothelial',['endo','vascular','lymphatic']),('Neural/glial',['neuro','neuron','glial','schwann']),('Myeloid',['myeloid','mac','mono','kupffer','dendritic','cdc','pdc','neutro','granulo','mast']),('Lymphoid',['lymph','t/nk','cd4','cd8',' nk','b/plasma','b cell','plasma']),('Epithelial/parenchymal',['appendage','epidermis','kc','keratin','epithelial','alveolar','at1','at2','basal','ciliated','hepato','cholang','tubule','duct','loop','collecting','cardiomyocyte','epicardial'])]
    return next((label for label,words in rules if any(w in text for w in words)),'Other')
d['display_label']=d.apply(reference,axis=1)
def raster(size,focus=None):
    xx=np.clip(np.rint((xy[:,0]-limits[0,0])/(limits[1,0]-limits[0,0])*(size-1)).astype(int),0,size-1)
    yy=np.clip(np.rint((limits[1,1]-xy[:,1])/(limits[1,1]-limits[0,1])*(size-1)).astype(int),0,size-1)
    im=np.full((size,size,3),255,np.uint8)
    mask=np.ones(len(d),bool) if focus is None else focus
    if focus is not None:im[yy,xx]=[235,235,235]
    for label,c in PALETTE.items():
        keep=mask&d.display_label.eq(label).to_numpy();im[yy[keep],xx[keep]]=[int(c[j:j+2],16) for j in [1,3,5]]
    return im
def show(ax,img):ax.imshow(img,interpolation='nearest',rasterized=True);ax.set_axis_off();ax.set_aspect('equal')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
st=pd.read_csv(PREV/'tables/SOX9_target_expression_statistics_exploratory.csv');st.to_csv(OUT/'tables/SOX9_statistics_unchanged.csv',index=False,encoding='utf-8-sig')
tabs={organ:pd.read_csv(PREV/'tables'/f'{organ}_target_sample_unit_expression.csv') for organ in ORDER}
fig=plt.figure(figsize=(210/25.4,75/25.4));W,H=210,75
fig.text(1/W,73/H,'b',fontsize=12,fontweight='bold',va='top');fig.text(84/W,73/H,'c',fontsize=12,fontweight='bold',va='top')
ax=fig.add_axes([2/W,9/H,57/W,57/H]);show(ax,raster(1700))
leg=fig.add_axes([60/W,27/H,24/W,38/H]);leg.set_axis_off();leg.set_xlim(0,1);leg.set_ylim(0,1)
for i,((label,col),name) in enumerate(zip(PALETTE.items(),['Epi/parenchymal','Fib/stromal','Endothelial','Neural/glial','Myeloid','Lymphoid','Other'])):
    y=.96-i*.15;leg.scatter([.02],[y],s=12,color=col,clip_on=False);leg.text(.105,y,name,fontsize=8,va='center')
# Same expressions, statistical units, P/FDR and frozen coordinates; this is a layout-only revision.
for i,organ in enumerate(ORDER):
    center=98+i*24;left=center-10
    fig.text(center/W,68/H,organ.title(),fontsize=9,weight='bold',ha='center')
    ax=fig.add_axes([left/W,42/H,20/W,20/H]);show(ax,raster(700,d.organ.eq(organ).to_numpy()))
    fig.text(center/W,35/H,NAMES[organ],fontsize=8,weight='bold',ha='center')
    ax=fig.add_axes([left/W,10/H,20/W,20/H]);row=st[st.organ==organ].iloc[0];tab=tabs[organ]
    v=[tab.loc[tab.condition_pair==g,'SOX9_mean_expression'].to_numpy() for g in [row.reference_group,row.case_group]]
    ax.bar([0,1],[x.mean() for x in v],width=.58,color=['#CBD6E2','#F2E8CA'],edgecolor='black',linewidth=.6)
    rng=np.random.default_rng(20261009+i)
    for j,x in enumerate(v):ax.scatter(j+rng.uniform(-.1,.1,len(x)),x,s=5,color='#333333',alpha=.7,linewidths=0,zorder=4)
    ax.set_ylim(0,3.05);ax.set_yticks([0,1,2,3]);ax.set_xticks([0,1],['Ref','Fib']);ax.tick_params(length=1.7,pad=1)
    if i==0:ax.set_ylabel(r'$SOX9$ expression',labelpad=2)
    else:ax.tick_params(labelleft=False)
    ax.text(.5,.965,f'P={row.P_exploratory:.3g}\nFDR={row.FDR_exploratory_BH_5:.3g}',ha='center',va='top',fontsize=8,transform=ax.transAxes)
fig.text(4/W,2.0/H,'P/FDR: exploratory sample-unit tests. Original UMAP coordinates retained.',fontsize=8)
fig.canvas.draw();sizes=[t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text()];assert min(sizes)>=8 and max(sizes)<=12
stem=OUT/'figures/Figure1bc_horizontal_selected_controls'
for ext in ['png','pdf','svg']:fig.savefig(stem.with_suffix('.'+ext),dpi=600,facecolor='white')
plt.close(fig)

# Table S1: only the 144 source-library records actually represented in the current plotted cell set.
records=[]
for _,r in m.iterrows():
    a=d[d.record_id==r.record_id];assert len(a)>0
    organ=r.organ;target=a[a.bulk_annotation=='Epithelial/parenchymal'] if organ=='liver' else a[a.cell_annotation.isin(TARGETS[organ])]
    n=len(target);total=float(target.SOX9_nonnegative_expression.sum())
    missing='Target cell population absent; mean is not evaluable' if n==0 else ''
    warnings=[]
    if pd.notna(r.warnings):warnings.append(r.warnings)
    if 0<n<20:warnings.append('Fewer than 20 target cells; imprecise expression estimate')
    group='; '.join(sorted(a.condition_pair.unique()))
    records.append({'Record ID':r.record_id,'Organ':organ,'Species':r.species,'Dataset':r.dataset,'Sample/library ID':r['sample'],'Patient ID as recorded':r.patient_id_source,'True tissue':r.actual_tissue,'Clinical condition':r.clinical_group,'Figure comparison group':group,'Timepoint':r.timepoint,'Age/stage':r.age_stage,'Sex':r.sex,'Anatomic site':r.anatomical_site,'Cell enrichment':r.cell_enrichment,'Source accession':r.source_accession,'Assay':'scRNA-seq','Figure b cells':len(a),'Target cell population':NAMES[organ],'Target cells':n,'Target SOX9 expression sum':total,'Target SOX9 mean expression':total/n if n else None,'Expression statistics inclusion':'Included' if n else 'Not evaluable','Missing target reason':missing,'Source Fib compartment cells':int(r.Fib_compartment_cells),'Source Fib SOX9 mean expression':r.Fib_compartment_SOX9_mean,'Plot unit(s)':'; '.join(sorted(a.sample_unit.unique())),'Expression basis':'; '.join(sorted(a.expression_basis.unique())),'Metadata mapping status':r.mapping_status,'Source panel':r.source_panel,'Warnings':'; '.join(warnings)})
s1=pd.DataFrame(records).sort_values('Record ID').reset_index(drop=True)
assert set(s1['Record ID'])==set(d.record_id) and s1['Figure b cells'].sum()==312292
assert (s1['Target cells']<=s1['Figure b cells']).all()
for organ in ORDER:
    a=s1[s1.Organ==organ];assert a['Target cells'].sum()==tabs[organ].target_cells.sum()
    previous=pd.read_csv(PREV/'tables/target_source_library_expression.csv');p=previous[previous.organ==organ].set_index('record_id')
    for _,r in a[a['Target cells']>0].iterrows():assert np.isclose(r['Target SOX9 mean expression'],p.loc[r['Record ID'],'SOX9_mean_expression'],atol=1e-12)
s1.to_csv(OUT/'tables/Table_S1_current_Figure1bc_samples.csv',index=False,encoding='utf-8-sig')
definitions=pd.DataFrame([
    ['Scope','Only source-library records actually used for current Figure1b-c. Excluded historical/scalp references are not listed.'],
    ['Observation unit','One source-library record. Library count is not independent patient/animal count.'],
    ['Included samples','144 records: Skin 67, Lung 43, Kidney 5, Heart 9, Liver 20; 312,292 plotted cells.'],
    ['Seven skin controls','F086, F087, F089, F090, F091, F092, F093. Source names preserved.'],
    ['Four controls without Fib_Myo','F086, F087, F092 and F093 remain in the UMAP and this table. Fib_Myo mean is unavailable, not zero.'],
    ['Target cells','Skin Fib_Myo; lung AT epithelium; kidney collecting duct epithelium; heart Fib/activated Fib; liver original Epi/Par broad group.'],
    ['Expression','Average nonnegative normalized RNA within target cells, including zeros. Not SOX9-positive fraction.'],
    ['Normalization','Skin source normalized X; public source SOX9 counts / original full-cell n_counts, log1p(CP10K). Values are not directly comparable across sources.'],
    ['Statistical units','Figure c points use legacy sample_unit means. The table gives source-library-level means. These aggregation levels differ.'],
    ['P and FDR','Unchanged exploratory two-sided Welch tests of legacy sample-unit means; BH FDR across five comparisons. Biological independence unverified.'],
    ['Scar-side membership','Kept unchanged, including acute wounds, scarless healing and five PBMC libraries with no Fib_Myo cells. True tissue and clinical condition are shown.'],
    ['Missing metadata','Blank source fields are displayed as NR (not recorded). No ages, sex or independent animal IDs inferred.'],
    ['UMAP','Original integrated coordinates retained. No new integration, annotation or clustering in this layout revision.'],
],columns=['Definition','Details'])
def payload(name,f):return {'name':name,'headers':list(f.columns),'rows':json.loads(f.replace([np.inf,-np.inf],np.nan).to_json(orient='values',force_ascii=False))}
# NR applies to source metadata fields only. Preserve numerical zero, unavailable target means and genuine empty warning fields.
for c in ['Patient ID as recorded','Timepoint','Age/stage','Sex','Anatomic site','Cell enrichment','Source accession','Metadata mapping status']:
    s1[c]=s1[c].fillna('NR').replace('', 'NR')
(OUT/'audit/workbook_payload.json').write_text(json.dumps([payload('Table S1',s1),payload('Definitions',definitions)],ensure_ascii=False,indent=2),encoding='utf8')
qa={'current_sample_records':len(s1),'current_cells':len(d),'organ_records':s1.Organ.value_counts().to_dict(),'expression_statistics_library_records':int((s1['Target cells']>0).sum()),'missing_target_records':s1.loc[s1['Target cells']==0,'Record ID'].tolist(),'metadata_scope':'Current figure only; no historical excluded records','layout_mm':[210,75],'font_family':'Arial','font_range_pt':[8,12],'analysis_membership_and_statistics_unchanged':True,'skin_target_remains_Fib_Myo':True,'zero_imputation_of_missing_target':False}
(OUT/'audit/scope_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Figure1 b-c横排及Table S1</title><style>body{{font:15px Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#193b52;background:#f4f7f9}}section{{background:white;padding:22px;margin:16px 0}}img{{max-width:100%;height:auto}}table{{border-collapse:collapse;white-space:nowrap;font-size:12px}}td,th{{padding:8px;border-bottom:1px solid #ddd}}th{{background:#193e59;color:white;position:sticky;top:0}}.scroll{{overflow:auto;max-height:700px}}.warn{{border-left:5px solid #c45b48}}</style><h1>Figure 1 b–c 横排与当前数据 Table S1</h1><section><img src="../figures/Figure1bc_horizontal_selected_controls.png"><p><a href="../figures/Figure1bc_horizontal_selected_controls.pdf">PDF</a> · <a href="../figures/Figure1bc_horizontal_selected_controls.svg">SVG</a> · <a href="../outputs/Table_S1_current_Figure1bc_samples.xlsx">Table S1 Excel</a> · <a href="../tables/Table_S1_current_Figure1bc_samples.csv">Table S1 CSV</a></p></section><section class="warn"><p>本次仅调整版式并整理metadata，不改变上一版样本、对应细胞群或统计结果。皮肤仍为Fib_Myo。7个对照全部保留于UMAP及S1，其中4个没有Fib_Myo，表达均值记为不可计算，不填零。尚未确认改为全Fib口径。</p><p>S1只列当前作图的144条来源文库记录，共312,292个细胞，不含历史排除数据。5条PBMC文库因“旧瘢痕侧不变”仍在本图中，S1保留其真实组织，不能视为皮肤独立样本。急性创面/无瘢痕及成年/P8等混合也保留。P/FDR为探索性旧统计单位结果，不作为患者层面正式推断。</p><p>b采用颜色和图例注释，未复制参考图未经核实位置的文字标签。c柱为SOX9表达均值而非阳性比例；小图共享0–3.05纵轴。Arial 8–12pt，整幅210×75mm，PNG 600dpi及可编辑SVG/PDF文字。</p></section><section><h2>Table S1 当前使用样本</h2><p>筛选：<input id="search" placeholder="器官、样本、组织或临床组"></p><div class="scroll" id="sampletable">{s1.to_html(index=False,na_rep='NA',float_format=lambda x:f'{x:.4f}',escape=True)}</div></section><section><h2>定义</h2>{definitions.to_html(index=False,escape=True)}</section><script>document.getElementById('search').addEventListener('input',e=>{{document.querySelectorAll('#sampletable tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(e.target.value.toLowerCase())?'':'none')}});</script></html>'''
(OUT/'report/Figure1bc_horizontal_TableS1.html').write_text(page,encoding='utf8')
print(json.dumps(qa,ensure_ascii=False,indent=2))
