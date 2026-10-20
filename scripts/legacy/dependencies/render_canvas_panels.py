from pathlib import Path
import json, sys, ast, hashlib, importlib.util, io
import h5py, numpy as np, pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PIL import Image
import fitz
sys.stdout.reconfigure(encoding='utf8')
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path('data/external/INPUT_6ab40abaef/SOX9投稿20261008')
OLD=Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
SKILL=Path('data/external/INPUT_8c05c810a6/run_figure2_ae_universal.py')
SPATIAL_REF=Path('data/external/INPUT_214c33b60e/pasted-text.txt')
RATROOT=Path('data/external/INPUT_01762fefba/spatial_h5ad_compatible')
ROIFILE=Path('data/external/INPUT_d28991e41f/FIGURE_INDEX.tsv')
SLOTS={2:[.006,.008,.752,.256],3:[.007,.719,.476,.962]}
PAGEPIX=(4961,7016)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def readcol(g):
    if isinstance(g,h5py.Dataset):return g.asstr()[:] if g.dtype.kind in ['O','S','U'] else g[:]
    if 'values' in g:return readcol(g['values'])
    c=g['codes'][:];v=readcol(g['categories']);a=v[np.maximum(c,0)].astype(object);a[c<0]='NA';return a
def rectangle(figure):
    with Image.open(SOURCE/f'Figure {figure}.tif') as im:
        native=im.convert('RGB');scaled_h=round(native.height/native.width*PAGEPIX[0]);native=native.resize((PAGEPIX[0],scaled_h),Image.Resampling.LANCZOS)
    base=Image.new('RGB',PAGEPIX,'white');base.paste(native,(0,0))
    f=SLOTS[figure];box=(round(f[0]*PAGEPIX[0]),round(f[1]*scaled_h),round(f[2]*PAGEPIX[0]),round(f[3]*scaled_h))
    W=(box[2]-box[0])*210/PAGEPIX[0];H=(box[3]-box[1])*297/PAGEPIX[1]
    return base,box,W,H
def style():
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
def export(fig,name):
    fig.canvas.draw();sizes=[x.get_fontsize() for x in fig.findobj(matplotlib.text.Text) if x.get_text()];assert min(sizes)>=8 and max(sizes)<=12
    for ext in ['pdf','svg','png']:fig.savefig(ROOT/'figures'/f'{name}.{ext}',dpi=600,facecolor='white')
    with Image.open(ROOT/'figures'/f'{name}.png') as im:im.convert('RGB').save(ROOT/'figures'/f'{name}.tiff',dpi=(600,600),compression='tiff_lzw')
    plt.close(fig)

# Use the immutable skill functions and its exact marker-dotplot body. No analytical/rendering formula is substituted.
assert sha(SKILL)==(ROOT/'audit/skill_renderer_sha256.txt').read_text().strip()
spec=importlib.util.spec_from_file_location('skill_figure2',SKILL);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
source_text=SKILL.read_text(encoding='utf8');tree=ast.parse(source_text);main=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='main')
start=next(i for i,x in enumerate(main.body) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='dot' for t in x.targets))
end=next(i for i,x in enumerate(main.body[start:],start) if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='panel_contract' for t in x.targets))
dot_nodes=main.body[start:end]
snippet='\n'.join(ast.get_source_segment(source_text,x) or '' for x in dot_nodes)
(ROOT/'code/skill_dotplot_body_exact.txt').write_text(snippet,encoding='utf8')
run=ROOT/'analysis/human_universal';manifest=json.loads((run/'00_manifest/Figure2_A-E_manifest.json').read_text(encoding='utf8'))
geometry=pd.read_csv(Path(manifest['geometry_path']),sep='\t',float_precision='round_trip');coordinates=geometry[['UMAP1','UMAP2']].to_numpy();labels=geometry.cluster.to_numpy();clusters=manifest['clusters']
actual=hashlib.sha256(geometry[['barcode','UMAP1','UMAP2']].to_csv(sep='\t',index=False).encode()).hexdigest();assert actual==manifest['geometry_sha256']
marker=pd.read_csv(run/'05_tables/Figure2_marker_allowlist_summary.tsv',sep='\t',float_precision='round_trip');ledger=pd.read_csv(run/'05_tables/Figure2_selected_marker_ledger.tsv',sep='\t',float_precision='round_trip')
mandatory=['SOX9','PDGFRA','VIM','ACTA2'];displayed=list(dict.fromkeys(mandatory+ledger.gene.tolist()));assert (ledger.groupby('label').size()==3).all()
palette=dict(zip(clusters,plt.cm.tab20(np.linspace(0,1,len(clusters)))))
style();base,box,W,H=rectangle(2);fig=plt.figure(figsize=(W/25.4,H/25.4))
fig.text(1/W,(H-1)/H,'a',fontsize=12,weight='bold',va='top');fig.text(61/W,(H-1)/H,'b',fontsize=12,weight='bold',va='top')
fig.text(8/W,(H-1)/H,'Human Fib clusters',fontsize=9,weight='bold',va='top')
fig.text(68/W,(H-1)/H,'Marker expression',fontsize=9,weight='bold',va='top')
ax=fig.add_axes([3/W,19/H,48/W,37/H]);module.plot_cluster_umap(ax,coordinates,labels,palette,'');ax.set_xlabel('');ax.set_ylabel('');ax.set_aspect('equal')
highest=manifest['highest_target_cluster'];xy=coordinates[labels==highest];center=np.median(xy,axis=0)
ax.annotate(highest+'*',center,xytext=(8,7),textcoords='offset points',fontsize=9,weight='bold',color='black',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':.5},arrowprops={'arrowstyle':'-','linewidth':.6})
ax.legend(loc='upper left',bbox_to_anchor=(-.07,-.02),ncol=4,frameon=False,fontsize=8,markerscale=3,handletextpad=.2,columnspacing=.35,borderaxespad=0)
fig.text(4/W,1.2/H,highest+r'*: highest mean $SOX9$',fontsize=8,va='bottom')
dot_axis=fig.add_axes([77/W,11/H,67/W,40/H]);env={'fig':fig,'dot_axis':dot_axis,'dot_letter':'b','marker':marker,'displayed_genes':displayed,'cluster_order':clusters,'mandatory':mandatory,'np':np}
exec(compile(ast.Module(body=dot_nodes,type_ignores=[]),str(SKILL),'exec'),env)
dot_axis.set_title('',loc='left');dot_axis.set_xlabel('')
cax=fig.add_axes([78/W,55/H,27/W,1.8/H]);cb=fig.colorbar(ScalarMappable(norm=Normalize(-2,2),cmap='RdBu_r'),cax=cax,orientation='horizontal');cb.set_ticks([-2,0,2]);cb.ax.tick_params(labelsize=8,length=1,pad=.5)
fig.text(109/W,56/H,'Mean expression (z)',fontsize=8,va='center')
fig.text(77/W,1/H,'Dot size: percent expressing',fontsize=8,va='bottom')
export(fig,'Figure2AB_human_Fib_slot')
marker[marker.gene.isin(['SOX9','ACTA2','KRT14','COL1A1'])].to_csv(ROOT/'tables/Figure2_highest_cluster_identity_audit.tsv',sep='\t',index=False)

# Compile only the relevant functions from the user-supplied spatial renderer; skip all top-level filesystem and main actions.
assert sha(SPATIAL_REF)=='ef5037fc66836453b7b518a9b567908b24979adcc854ab0b42d768fc44c48d2c'
reference=ast.parse(SPATIAL_REF.read_text(encoding='utf8'));fnames={'rectangular_limits','format_square','format_zoom_out','add_zoom_roi_box','plot_distribution','nice_scale_bar','add_scale_bar'}
functions=[x for x in reference.body if isinstance(x,ast.FunctionDef) and x.name in fnames];assert len(functions)==7
space={'np':np,'Rectangle':Rectangle};exec(compile(ast.Module(body=functions,type_ignores=[]),str(SPATIAL_REF),'exec'),space)
rois=pd.read_csv(ROIFILE,sep='\t').set_index('sample')
sections=[('Adult','1 dpw','Adult_1d_C01830C1D1'),('Adult','5 dpw','Adult_5d1_C01830E6F6'),('Adult','19 dpw','Adult_19d_B03424F4'),('E16.5','6 h','E16.5_6H_C02938E1'),('E16.5','3 dpw','E16.5_3D_D03056D5'),('E16.5','6 dpw','E16.5_6D_C02938F4')]
oldmanifest=json.loads((OLD/'00_manifest/Figure3G_rat_spatial_manifest_v05.json').read_text(encoding='utf8'));expected={x['sample']:x['sha256'] for x in oldmanifest['inputs']}
base3,box3,W3,H3=rectangle(3);fig=plt.figure(figsize=(W3/25.4,H3/25.4));fig.text(.5/W3,(H3-1.5)/H3,'g',fontsize=12,weight='bold',va='top')
fig.text(5/W3,(H3-1)/H3,r'Fib: $Sox9>0$ and $Col1a1>0$',fontsize=8,va='top')
fig.text(7/W3,(H3-10)/H3,'Adult',fontsize=9,weight='bold');fig.text((W3/2+7)/W3,(H3-10)/H3,'E16.5',fontsize=9,weight='bold')
legend=fig.add_axes([54/W3,(H3-5)/H3,43/W3,4/H3]);legend.set_axis_off();legend.set_xlim(0,1);legend.set_ylim(0,1)
legend.scatter([.02],[.5],s=12,color='#FF1A1A');legend.text(.07,.5,'Target',fontsize=8,va='center');legend.scatter([.46],[.5],s=12,color='#E0E0E0');legend.text(.51,.5,'Other bins',fontsize=8,va='center')
inputs=[];tiles=[]
for i,(stage,time,sample) in enumerate(sections):
    path=RATROOT/f'{sample}_Annotated_Spatial.h5ad';before=(path.stat().st_size,path.stat().st_mtime_ns);digest=sha(path);assert digest==expected[sample],sample+' input changed'
    with h5py.File(path,'r') as h:
        vg=h['var'];names=readcol(vg[vg.attrs['_index']]);j=[np.flatnonzero(names==gene) for gene in ['Sox9','Col1a1']];assert all(len(v)==1 for v in j)
        sub=readcol(h['obs/sub_labels']).astype(str);fib=np.char.startswith(sub,'Fib_')
        key='X_rot_spatial' if 'X_rot_spatial' in h['obsm'] else 'spatial';xy=h['obsm'][key][:].astype(float)
        g=h['X'];assert g.attrs['encoding-type']=='csr_matrix';matrix=sparse.csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g.attrs['shape']))
        vals=matrix[:,[int(x[0]) for x in j]].toarray();assert np.isfinite(vals).all() and (vals>=0).all();target=fib&(vals>0).all(axis=1)
    r=rois.loc[sample];assert key==r.coordinate_key;roi=tuple(float(r[k]) for k in ['roi_xmin','roi_xmax','roi_ymin','roi_ymax'])
    col=i//3;row=i%3;x0=col*(W3/2)+5;rowheight=(H3-12)/3;bottom=H3-11-(row+1)*rowheight+.6
    fig.text((x0-3)/W3,(bottom+6)/H3,time,fontsize=8,rotation=90,va='center',ha='center')
    original_full=space['rectangular_limits'](xy)
    full=(min(original_full[0],roi[0]),max(original_full[1],roi[1]),min(original_full[2],roi[2]),max(original_full[3],roi[3]))
    am=fig.add_axes([x0/W3,bottom/H3,27/W3,12/H3]);az=fig.add_axes([(x0+29)/W3,bottom/H3,12/W3,12/H3])
    space['plot_distribution'](am,xy,target,full,'#FF1A1A',zoom=False);space['add_zoom_roi_box'](am,roi);bar=space['add_scale_bar'](am,full,.5)
    space['plot_distribution'](az,xy,target,roi,'#FF1A1A',zoom=True);zbar=space['add_scale_bar'](az,roi,.5)
    # Equivalent shorter units match the user's manuscript scale-bar labels, without changing bar lengths.
    if bar==1000:am.texts[-1].set_text('1 mm')
    if zbar==1000:az.texts[-1].set_text('1 mm')
    inroi=(xy[:,0]>=roi[0])&(xy[:,0]<=roi[1])&(xy[:,1]>=roi[2])&(xy[:,1]<=roi[3])
    rec={'stage':stage,'wound_time':time,'sample':sample,'source':str(path),'sha256':digest,'spatial_observations':len(sub),'source_Fib_observations':int(fib.sum()),'Sox9_positive_Fib':int(np.sum(fib&(vals[:,0]>0))),'Col1a1_positive_Fib':int(np.sum(fib&(vals[:,1]>0))),'Sox9_Col1a1_double_positive':int(target.sum()),'Sox9_only_Fib':int(np.sum(fib&(vals[:,0]>0)&(vals[:,1]==0))),'double_positive_cutoff1':int(np.sum(fib&(vals>1).all(axis=1))),'target_in_registered_zoom':int(np.sum(target&inroi)),'coordinate_key':key,'additional_rotation':0,'um_per_coordinate_unit':.5,'full_scale_bar_um':bar,'zoom_scale_bar_um':zbar,'roi_xmin':roi[0],'roi_xmax':roi[1],'roi_ymin':roi[2],'roi_ymax':roi[3],'independent_animal_id':'Not recorded','statistical_status':'One registered section per displayed time; descriptive only'}
    inputs.append(rec);assert before==(path.stat().st_size,path.stat().st_mtime_ns)
    print('SPATIAL',sample,'target',int(target.sum()),flush=True)
export(fig,'Figure3G_rat_spatial_slot');pd.DataFrame(inputs).to_csv(ROOT/'tables/Figure3G_used_sections_metadata_counts.csv',index=False,encoding='utf-8-sig')

# Preserve every experimental/reference pixel outside the explicitly registered panel slot.
composites=[]
for number,name,bg,bounds in [(2,'Figure2AB_human_Fib_slot',base,box),(3,'Figure3G_rat_spatial_slot',base3,box3)]:
    source_original=bg.copy();blank=Image.new('RGB',(bounds[2]-bounds[0],bounds[3]-bounds[1]),'white');bg.paste(blank,(bounds[0],bounds[1]))
    paneldoc=fitz.open(ROOT/'figures'/f'{name}.pdf');pdf=fitz.open();page=pdf.new_page(width=210/25.4*72,height=297/25.4*72)
    buffer=io.BytesIO();bg.save(buffer,format='PNG');page.insert_image(page.rect,stream=buffer.getvalue())
    rect=fitz.Rect(bounds[0]/PAGEPIX[0]*page.rect.width,bounds[1]/PAGEPIX[1]*page.rect.height,bounds[2]/PAGEPIX[0]*page.rect.width,bounds[3]/PAGEPIX[1]*page.rect.height)
    page.show_pdf_page(rect,paneldoc,0,keep_proportion=False)
    dest=ROOT/'figures'/f'Figure{number}_completed_A4';pdf.save(dest.with_suffix('.pdf'));pdf.close()
    pix=paneldoc[0].get_pixmap(matrix=fitz.Matrix((bounds[2]-bounds[0])/paneldoc[0].rect.width,(bounds[3]-bounds[1])/paneldoc[0].rect.height),alpha=False)
    panel=Image.frombytes('RGB',(pix.width,pix.height),pix.samples).resize(blank.size,Image.Resampling.LANCZOS);proof=source_original.copy();proof.paste(panel,(bounds[0],bounds[1]))
    original=np.asarray(source_original);current=np.asarray(proof);mask=np.ones(PAGEPIX[::-1],bool);mask[bounds[1]:bounds[3],bounds[0]:bounds[2]]=False
    assert np.array_equal(original[mask],current[mask])
    proof.save(dest.with_suffix('.png'),dpi=(600,600));proof.save(dest.with_suffix('.tiff'),dpi=(600,600),compression='tiff_lzw')
    preview=proof.copy();preview.thumbnail((1800,2500));preview.save(ROOT/'audit'/f'Figure{number}_completed_preview.png')
    composites.append({'figure':number,'slot_fraction':SLOTS[number],'slot_pixels_in_A4':list(bounds),'slot_mm':[(bounds[2]-bounds[0])*210/PAGEPIX[0],(bounds[3]-bounds[1])*297/PAGEPIX[1]],'outside_slot_pixel_changes':0,'experimental_panels':'Inherited unchanged from user TIFF; raw experimental measurements and font sizes not independently verified'})
qa={'human_fresh_skill_execution':True,'skill_source_sha256':sha(SKILL),'geometry_sha256':actual,'highest_SOX9_cluster':highest,'n_Fib_cells':len(geometry),'n_clusters':len(clusters),'C8_ACTA2_mean':float(marker[(marker.cluster==highest)&(marker.gene=='ACTA2')].mean_expression.iloc[0]),'C8_KRT14_mean':float(marker[(marker.cluster==highest)&(marker.gene=='KRT14')].mean_expression.iloc[0]),'myofibroblast_identity_not_assumed':True,'spatial_marker_gate':'Source Fib_ annotation AND Sox9>0 AND Col1a1>0','spatial_sections':len(inputs),'spatial_renderer_functions_unchanged':True,'reference_spatial_source_sha256':sha(SPATIAL_REF),'scale_provenance':'0.5 um/unit inherited from user reference renderer and original rat 30-section manifest','zoom_ROI_selection':'Registered original ROIs reused without enrichment retuning for new markers','new_spatial_group_test':'Not run; n=1 section per displayed stage-time, no independent animal mapping','composites':composites}
(ROOT/'audit/panel_execution_and_canvas.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8');print('CANVAS PANELS AND A4 PROOFS COMPLETE',flush=True)
