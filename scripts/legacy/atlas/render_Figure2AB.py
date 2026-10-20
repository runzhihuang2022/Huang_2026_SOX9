from pathlib import Path
import ast, importlib.util, hashlib, json, shutil, sys, io
import numpy as np, pandas as pd, h5py, anndata as ad
from scipy import sparse
import fitz
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_full_selected import ROOT, PROJECT, SKILL, sha
sys.stdout.reconfigure(encoding='utf8')

def main():
    original=PROJECT/'20261009_Figure2_3_canvas_completion_v01/code/render_canvas_panels.py'
    source_text=original.read_text(encoding='utf8')
    prefix=source_text.split('# Compile only the relevant functions from the user-supplied spatial renderer',1)[0]
    assert source_text.count('# Compile only the relevant functions from the user-supplied spatial renderer')==1
    namespace={'__file__':str(__file__),'__name__':'selected19_Figure2_canvas'}
    # Execute the approved A/B canvas and pristine skill functions with the current versioned input root only.
    exec(compile(prefix,str(original),'exec'),namespace)
    fig=namespace['fig'];ax=namespace['ax'];dot_axis=namespace['dot_axis'];W,H=namespace['W'],namespace['H']
    ax.set_position([3/W,23/H,48/W,33/H])
    ax.legend(loc='upper left',bbox_to_anchor=(-.07,-.02),ncol=5,frameon=False,fontsize=8,markerscale=3,handlelength=.8,handletextpad=.25,columnspacing=.55,borderaxespad=0)
    for tick in dot_axis.get_xticklabels():tick.set_rotation(90);tick.set_horizontalalignment('center')
    namespace['export'](fig,'Figure2AB_human_Fib_slot')
    shutil.copy2(original,ROOT/'code/Figure2_reference_canvas_PRISTINE.py')
    (ROOT/'audit/Figure2_allowed_layout_changes.json').write_text(json.dumps({'source':str(original),'source_sha256':sha(original),'allowed_input_substitution':'Current versioned ROOT; full20,270 Fib input and geometry','allowed_layout_reflow':['17-cluster legend: 5 columns with compact handles; UMAP lower edge lifted 4mm to separate four legend rows from footer','x cluster labels: 90 degrees'],'analytical_functions_and_exact_skill_dotplot_body_unchanged':True,'KRT14_removed_by_marker_allowlist_not_by_expression_deletion':True,'canvas_mm':[W,H]},ensure_ascii=False,indent=2),encoding='utf8')
    geometry=namespace['geometry'];labels=geometry.cluster.to_numpy();highest=namespace['highest'];coords=namespace['coordinates'];palette=namespace['palette'];module=namespace['module']
    g=geometry.set_index('barcode')
    source=ROOT/'01_input_audit/Selected19_full_QC_Fib_normalized_RNA.h5ad'
    a=ad.read_h5ad(source)
    assert set(a.obs_names)==set(g.index)
    g=g.reindex(a.obs_names);cluster=g.cluster.to_numpy();gene_list=['SOX9','COL1A1','COL1A2','COL3A1','ACTA2','VIM','PDGFRA','KRT14','KRT19','EPCAM','PTPRC','PECAM1','VWF','RGS5','CSPG4']
    rows=[]
    for c in namespace['clusters']:
        mask=cluster==c
        for gene in gene_list:
            if gene not in a.var_names:continue
            vals=a.X[mask,a.var_names.get_loc(gene)].toarray().ravel()
            rows.append({'cluster':c,'gene':gene,'n_cells':int(mask.sum()),'mean_expression':float(vals.mean()),'percent_expressed':float((vals>0).mean()*100),'displayed_in_main_dotplot':gene in namespace['displayed']})
    purity=pd.DataFrame(rows);purity.to_csv(ROOT/'tables/Figure2_all_cluster_marker_and_purity_audit.csv',index=False,encoding='utf-8-sig')
    sox=purity[purity.gene.eq('SOX9')].sort_values(['mean_expression','percent_expressed'],ascending=False)
    sox.to_csv(ROOT/'tables/Figure2_cluster_SOΧ9_ranked.csv',index=False,encoding='utf-8-sig')
    for cols,name in [(['sub_labels'],'Figure2_cluster_original_subtype_counts.csv'),(['record_id','clinical_group'],'Figure2_cluster_library_counts.csv')]:
        obs=a.obs.copy();obs['cluster']=cluster;obs.groupby(['cluster']+cols,observed=True).size().reset_index(name='cells').to_csv(ROOT/'tables'/name,index=False,encoding='utf-8-sig')
    hi=purity[purity.cluster.eq(highest)]
    kr=hi[hi.gene.eq('KRT14')];pg=hi[hi.gene.eq('PDGFRA')]
    boundary='Highest observed SOX9 within the supplied original Fibroblast gate; unsupervised label only, no definitive myofibroblast identity.'
    if len(kr) and float(kr.percent_expressed.iloc[0])>25:boundary+=' KRT14-high cells require epithelial contamination/doublet/state review; removing KRT14 from display does not resolve this uncertainty.'
    qa={'status':'RENDERED_REQUIRES_FINAL_VISUAL_QA','full_Fib_cells':len(a),'skin_libraries':int(a.obs.record_id.nunique()),'n_clusters':len(namespace['clusters']),'highest_SOΧ9_cluster':highest,'highest_cluster_cells':int((cluster==highest).sum()),'highest_markers':hi.to_dict('records'),'highest_source_subtypes':a.obs.loc[cluster==highest,'sub_labels'].astype(str).value_counts().to_dict(),'highest_library_counts':a.obs.loc[cluster==highest,'record_id'].value_counts().to_dict(),'biological_boundary':boundary,'geometry_sha256':namespace['actual'],'shared_Figure1_skin_library_ids':sorted(a.obs.record_id.unique()),'new_panel_fonts_pt':[8,12],'main_panel_canvas_mm':[W,H]}
    (ROOT/'audit/Figure2_CLUSTER_IDENTITY_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    # Auxiliary target-feature and highest-cluster views on the exact frozen skill UMAP.
    namespace['style']()
    fig,axes=plt.subplots(1,2,figsize=(210/25.4,90/25.4));fig.subplots_adjust(left=.045,right=.97,top=.86,bottom=.19,wspace=.18)
    module.plot_cluster_umap(axes[0],coords,labels,palette,'Highest SOX9: '+highest,highlight=highest)
    target=a.X[:,a.var_names.get_loc('SOX9')].toarray().ravel();order=np.argsort(target)
    # Geometry input is reindexed to the RNA order; use matching coordinate rows for feature values.
    cg=g[['UMAP1','UMAP2']].to_numpy()
    im=axes[1].scatter(cg[order,0],cg[order,1],c=target[order],s=1.4,cmap='viridis',linewidths=0,rasterized=True)
    module.clean_umap(axes[1]);axes[1].set_title(r'$SOX9$ expression',loc='left',fontsize=10)
    cb=fig.colorbar(im,ax=axes[1],fraction=.035,pad=.02);cb.ax.tick_params(labelsize=8);cb.set_label('Source normalized RNA',fontsize=8)
    fig.text(.03,.025,'C labels are unsupervised. Highest observed expression is not proof of a specific biological subtype.',fontsize=8)
    namespace['export'](fig,'Figure2_SOΧ9_highest_and_feature_shared_geometry')
    # Rebuild the manuscript A4 proof by replacing only the registered A/B slot.
    before=namespace['base'].copy();bounds=namespace['box'];blank=before.copy();blank.paste(Image.new('RGB',(bounds[2]-bounds[0],bounds[3]-bounds[1]),'white'),bounds[:2])
    panel=fitz.open(ROOT/'figures/Figure2AB_human_Fib_slot.pdf');pdf=fitz.open();page=pdf.new_page(width=210/25.4*72,height=297/25.4*72)
    buffer=io.BytesIO();blank.save(buffer,format='PNG');page.insert_image(page.rect,stream=buffer.getvalue())
    ww,hh=before.size;rect=fitz.Rect(bounds[0]/ww*page.rect.width,bounds[1]/hh*page.rect.height,bounds[2]/ww*page.rect.width,bounds[3]/hh*page.rect.height)
    page.show_pdf_page(rect,panel,0,keep_proportion=False)
    pdf.save(ROOT/'figures/Figure2_complete_A4_selected19_full_Fib.pdf',garbage=4,deflate=True)
    direct=before.copy();slot=Image.open(ROOT/'figures/Figure2AB_human_Fib_slot.png').convert('RGB').resize((bounds[2]-bounds[0],bounds[3]-bounds[1]),Image.Resampling.LANCZOS);direct.paste(slot,bounds[:2])
    mask=np.ones((hh,ww),dtype=bool);mask[bounds[1]:bounds[3],bounds[0]:bounds[2]]=False
    assert np.array_equal(np.asarray(before)[mask],np.asarray(direct)[mask])
    direct.save(ROOT/'figures/Figure2_complete_A4_selected19_full_Fib.png',dpi=(600,600));direct.save(ROOT/'figures/Figure2_complete_A4_selected19_full_Fib.tiff',dpi=(600,600),compression='tiff_lzw')
    (ROOT/'audit/Figure2_source_slot_preservation.json').write_text(json.dumps({'slot_pixels':bounds,'slot_mm':[W,H],'outside_slot_changed_pixels':0,'unchanged_experimental_content_source':str(namespace['SOURCE']/'Figure 2.tif'),'experimental_panel_measurements_not_recomputed':True},indent=2),encoding='utf8')
    print(json.dumps(qa,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
