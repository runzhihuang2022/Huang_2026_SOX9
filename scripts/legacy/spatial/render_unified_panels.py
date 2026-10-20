from __future__ import annotations

import ast
import hashlib
import html
import io
import json
import sys
from pathlib import Path

import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf8')
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path('data/external/INPUT_6ab40abaef/SOX9投稿20261008')
CURRENT3=Path('data/external/INPUT_743ea21501/Figure3_completed_A4.png')
BC=Path('data/external/INPUT_72ca0088ee/Figure1bc_horizontal_selected_controls.png')
REFERENCE=Path('data/external/INPUT_214c33b60e/pasted-text.txt')
PAGEPIX=(4961,7016)
SLOTS={1:[.642,.354,.985,.566],3:[.007,.719,.476,.962]}
YELLOW='#E6B800'
RED='#FF1A1A'
GREEN='#00A83B'
GRAY='#E0E0E0'
MANIFEST=json.loads((ROOT/'audit/INPUT_AND_MARKER_AUDIT.json').read_text(encoding='utf8'))
OUTPUTS=[]; BARS=[]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def load_functions():
    assert sha(REFERENCE)=='ef5037fc66836453b7b518a9b567908b24979adcc854ab0b42d768fc44c48d2c'
    source=REFERENCE.read_text(encoding='utf8'); tree=ast.parse(source)
    names={'rectangular_limits','format_square','format_zoom_out','add_zoom_roi_box','plot_distribution','nice_scale_bar','add_scale_bar'}
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==7
    space={'np':np,'Rectangle':Rectangle};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(REFERENCE),'exec'),space)
    (ROOT/'code/reference_spatial_functions_pristine.py').write_text('\n\n'.join(ast.get_source_segment(source,n) for n in nodes),encoding='utf8')
    return space


REF=load_functions()
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})


def data(name):
    with np.load(ROOT/'data'/f'{name}.npz') as n:
        return {k:n[k] for k in n.files}


def geometry(number):
    with Image.open(SOURCE/f'Figure {number}.tif') as im:
        sh=round(im.height/im.width*PAGEPIX[0])
        if number==1:
            im=im.convert('RGB').resize((PAGEPIX[0],sh),Image.Resampling.LANCZOS)
            base=Image.new('RGB',PAGEPIX,'white');base.paste(im,(0,0))
        else:
            with Image.open(CURRENT3) as current:base=current.convert('RGB')
            assert base.size==PAGEPIX
    s=SLOTS[number]; box=(round(s[0]*PAGEPIX[0]),round(s[1]*sh),round(s[2]*PAGEPIX[0]),round(s[3]*sh))
    W=(box[2]-box[0])*210/PAGEPIX[0];H=(box[3]-box[1])*297/PAGEPIX[1]
    return base,box,W,H


def save(fig,name):
    fig.canvas.draw();texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_text()]
    sizes=[t.get_fontsize() for t in texts]
    assert min(sizes)>=8 and max(sizes)<=12,(name,min(sizes),max(sizes))
    # Clip check focuses on labels and legends, not on deliberately cropped scatter points.
    renderer=fig.canvas.get_renderer(); bounds=fig.bbox
    clipped=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if b.x0<bounds.x0-1 or b.y0<bounds.y0-1 or b.x1>bounds.x1+1 or b.y1>bounds.y1+1:
            clipped.append({'text':t.get_text(),'bbox':list(b.bounds)})
    assert not clipped,(name,clipped)
    for ext in ['pdf','svg','png']:
        fig.savefig(ROOT/'figures'/f'{name}.{ext}',dpi=600,facecolor='white')
    with Image.open(ROOT/'figures'/f'{name}.png') as im:
        im.convert('RGB').save(ROOT/'figures'/f'{name}.tiff',dpi=(600,600),compression='tiff_lzw')
    OUTPUTS.append({'name':name,'mm':(fig.get_size_inches()*25.4).tolist(),'fonts_pt':[min(sizes),max(sizes)],'font_family':'Arial','dpi':600,'out_of_canvas_text':clipped})
    plt.close(fig)


def add_legend(fig,W,H,x,y):
    ax=fig.add_axes([x/W,y/H,40/W,4/H]);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
    ax.scatter([.015],[.5],s=11,c=YELLOW,edgecolor='#9B7C00',linewidth=.3,clip_on=False)
    ax.text(.065,.5,'Dual positive',fontsize=8,va='center')
    ax.scatter([.54],[.5],s=11,c=GRAY,clip_on=False);ax.text(.59,.5,'Other bins',fontsize=8,va='center')


def marker_title(fig,W,H,letter,rat=False):
    fig.text(.5/W,(H-1)/H,letter,fontsize=12,weight='bold',va='top')
    fig.text(5/W,(H-1)/H,'RNA:',fontsize=8,va='top')
    fig.text(13/W,(H-1)/H,'Col1a1' if rat else 'COL1A1',fontstyle='italic',color=GREEN,fontsize=8,va='top')
    fig.text(25.5/W,(H-1)/H,'+',fontsize=8,va='top')
    fig.text(29/W,(H-1)/H,'Sox9' if rat else 'SOX9',fontstyle='italic',color=RED,fontsize=8,va='top')


def scale_bar(ax,limits):
    bar=REF['add_scale_bar'](ax,limits,.5)
    label=ax.texts[-1]
    # Keep the original physical bar length and 8 pt label; left-align the label
    # above its bar so narrow, manuscript-sized zoom panels do not spill left.
    label.set_ha('left')
    label.set_position((limits[0]+(limits[1]-limits[0])*.06,label.get_position()[1]))
    if bar==1000:label.set_text('1 mm')
    return bar


def maps(ax,az,name,d,roi,add_boundary=None):
    xy=d['xy'];target=d['fib']&(d['vals'][:,:2]>0).all(axis=1)
    full0=REF['rectangular_limits'](xy)
    full=(min(full0[0],roi[0]),max(full0[1],roi[1]),min(full0[2],roi[2]),max(full0[3],roi[3]))
    REF['plot_distribution'](ax,xy,target,full,YELLOW,zoom=False)
    REF['add_zoom_roi_box'](ax,roi)
    bar=scale_bar(ax,full)
    REF['plot_distribution'](az,xy,target,roi,YELLOW,zoom=True)
    zbar=scale_bar(az,roi)
    if add_boundary is not None:
        (x0,y0),(x1,y1)=np.asarray(add_boundary)
        # The saved anatomical line is drawn only across its recorded y support.
        ax.plot([x0,x1],[y0,y1],color='#555555',linestyle=':',lw=.7,zorder=7)
    BARS.append({'name':name,'full_scale_bar_um':bar,'zoom_scale_bar_um':zbar,'um_per_coordinate_unit':.5,'roi':list(roi),'coordinate_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'positive_mask_sha256':hashlib.sha256(target.tobytes()).hexdigest(),'double_positive':int(target.sum()),'roi_selection_changed':False})


def human_regions():
    rows=[]
    for rec in MANIFEST['human_source_parts']:
        name=rec['data_key'];d=data(name);target=d['fib']&(d['vals'][:,:2]>0).all(axis=1)
        section='6 mph' if name.startswith('6mph') else rec['section']
        for region in ['Scarless side','Scar side']:
            m=d['fib']&(d['region']==region);n=int(m.sum());k=int((m&target).sum())
            rows.append({'section':section,'source_part':name,'region':region,'n_source_Fib':n,'dual_positive':k,'percent':100*k/n if n else None,'unit':'Spatial observations / bins','P_value':'NE','BH_FDR':'NE','test_status':'No biological-replicate test: donor independence unavailable'})
    bypart=pd.DataFrame(rows);bypart.to_csv(ROOT/'tables/Human_marker_by_source_part_region.csv',index=False,encoding='utf-8-sig')
    pooled=bypart.groupby(['section','region'],sort=False)[['n_source_Fib','dual_positive']].sum().reset_index()
    pooled['percent']=pooled.dual_positive/pooled.n_source_Fib*100
    pooled['P_value']='NE';pooled['BH_FDR']='NE'
    pooled.to_csv(ROOT/'tables/Human_marker_by_physical_section_region.csv',index=False,encoding='utf-8-sig')
    return pooled


def human_main(regions):
    base,box,W,H=geometry(1);fig=plt.figure(figsize=(W/25.4,H/25.4));marker_title(fig,W,H,'e')
    add_legend(fig,W,H,5,H-9)
    for i,name in enumerate(['2mph','6mph']):
        bottom=H-28.5-i*20.2
        fig.text(4/W,(bottom+15)/H,'2 mph' if i==0 else '6 mph',fontsize=9,weight='bold')
        d=data('main_'+name);roi=d['roi']
        am=fig.add_axes([3/W,bottom/H,23/W,13/H]);az=fig.add_axes([28/W,bottom/H,13/W,13/H])
        boundary=data('2mph')['boundary'] if i==0 else None
        maps(am,az,'human_main_'+name,d,roi,add_boundary=boundary)
    # Replace the old blank statistics placeholder with explicitly descriptive proportions.
    # Only the two displayed physical sections are shown. The 6 mph parts are pooled.
    ax=fig.add_axes([54/W,13/H,16/W,26/H])
    for i,section in enumerate(['2 mph','6 mph']):
        v=[regions.loc[(regions.section==section)&(regions.region==r),'percent'].iloc[0] for r in ['Scarless side','Scar side']]
        ax.plot([0,1],v,color='#888888',lw=.6,linestyle='-' if i==0 else '--')
        ax.scatter([0,1],v,marker='o' if i==0 else 's',s=11,c=YELLOW,edgecolor='#555555',linewidth=.4,zorder=4)
    ax.spines[['top','right']].set_visible(False);ax.set_ylim(0,regions[regions.section.isin(['2 mph','6 mph'])].percent.max()*1.5)
    ax.set_xticks([0,1],['SL','S']);ax.set_ylabel('Dual+ Fib (%)',labelpad=1);ax.tick_params(length=1.5,pad=1)
    ax.set_title('Descriptive',fontsize=8,pad=4)
    fig.text(53/W,4.5/H,'P=NE',fontsize=8);fig.text(53/W,1/H,'FDR=NE',fontsize=8)
    save(fig,'Figure1E_COLI_SOX9_unified_slot')
    return base,box


def rat_main():
    base,box,W,H=geometry(3);fig=plt.figure(figsize=(W/25.4,H/25.4));marker_title(fig,W,H,'g',rat=True)
    add_legend(fig,W,H,55,H-5)
    fig.text(7/W,(H-10)/H,'Adult',fontsize=9,weight='bold');fig.text((W/2+7)/W,(H-10)/H,'E16.5',fontsize=9,weight='bold')
    for i,rec in enumerate(MANIFEST['rat_main_sections']):
        name=rec['sample'];d=data(name);roi=d['roi'];col=i//3;row=i%3
        x0=col*(W/2)+5;rh=(H-12)/3;bottom=H-11-(row+1)*rh+.6
        fig.text((x0-3)/W,(bottom+6)/H,rec['wound_time'],fontsize=8,rotation=90,va='center',ha='center')
        am=fig.add_axes([x0/W,bottom/H,27/W,12/H]);az=fig.add_axes([(x0+29)/W,bottom/H,12/W,12/H])
        maps(am,az,name,d,roi)
    save(fig,'Figure3G_COLI_SOX9_unified_slot')
    return base,box


def channel_panel(ax,d,channel,color,roi=None):
    fib=d['fib'];xy=d['xy'];vals=d['vals'];target=fib&(vals[:,channel]>0)
    full=REF['rectangular_limits'](xy) if roi is None else roi
    REF['plot_distribution'](ax,xy,target,full,color,zoom=roi is not None)
    scale_bar(ax,full)


def channel_audits():
    # Real observed marker detection, not synthesized microscopy or protein intensity.
    for species,records in [('Human',MANIFEST['human_source_parts']),('Rat',MANIFEST['rat_main_sections'])]:
        W,H=210,297;fig=plt.figure(figsize=(W/25.4,H/25.4))
        fig.text(.03,.98,f'{species}: COL I / SOX9 marker concordance',fontsize=12,weight='bold',va='top')
        fig.text(.03,.95,'RNA detection in source Fib lineage; not immunofluorescence images.',fontsize=8,va='top')
        gene1='SOX9' if species=='Human' else 'Sox9';gene2='COL1A1' if species=='Human' else 'Col1a1'
        headers=[(38,gene1,RED),(105,gene2,GREEN),(174,'Dual positive',YELLOW)] if species=='Human' else [(27,gene1,RED),(76,gene2,GREEN),(125,'Dual positive',YELLOW),(174,'Registered zoom',YELLOW)]
        for x,text,c in headers:
            fig.text(x/W,.91,text,fontstyle='italic' if text in [gene1,gene2] else 'normal',color=c,fontsize=9,weight='normal' if text in [gene1,gene2] else 'bold',ha='center')
        rowheight=53 if species=='Human' else 34
        for i,rec in enumerate(records):
            name=rec['data_key'] if species=='Human' else rec['sample'];d=data(name)
            top=264-i*rowheight;bottom=top-rowheight+5
            label=rec['section'] if species=='Human' else rec['stage']+' '+rec['wound_time']
            fig.text(3/W,(top+1)/H,label,fontsize=8,weight='bold')
            mapheight=rowheight-9
            for j,(idx,c) in enumerate([(0,RED),(1,GREEN)]):
                left=5+j*69 if species=='Human' else 5+j*49
                width=62 if species=='Human' else 43
                ax=fig.add_axes([left/W,bottom/H,width/W,mapheight/H]);channel_panel(ax,d,idx,c)
            dual=d['fib']&(d['vals'][:,:2]>0).all(axis=1)
            ax=fig.add_axes([(143 if species=='Human' else 103)/W,bottom/H,(62 if species=='Human' else 43)/W,mapheight/H])
            REF['plot_distribution'](ax,d['xy'],dual,REF['rectangular_limits'](d['xy']),YELLOW)
            scale_bar(ax,REF['rectangular_limits'](d['xy']))
            if species=='Rat':
                roi=d['roi']
                REF['add_zoom_roi_box'](ax,roi)
                az=fig.add_axes([152/W,bottom/H,min(43,mapheight)/W,mapheight/H]);REF['plot_distribution'](az,d['xy'],dual,roi,YELLOW,zoom=True)
                scale_bar(az,roi)
        fig.text(.03,.06,'SOX9 channel: red. COL I-associated transcript channel: green. Dual detection: yellow.',fontsize=8)
        fig.text(.03,.045,'Gate: source Fib_ annotation and both RNA markers >0; grey = other spatial observations.',fontsize=8)
        fig.text(.03,.03,'Original microscope channels and measured pixel intensities are unchanged.',fontsize=8)
        save(fig,f'{species}_dual_marker_channel_audit_A4')


def compose(number,base,box,name):
    original=np.asarray(base);background=base.copy();background.paste('white',box)
    panel=fitz.open(ROOT/'figures'/f'{name}.pdf');doc=fitz.open();page=doc.new_page(width=210/25.4*72,height=297/25.4*72)
    b=io.BytesIO();background.save(b,format='PNG');page.insert_image(page.rect,stream=b.getvalue())
    rect=fitz.Rect(box[0]/PAGEPIX[0]*page.rect.width,box[1]/PAGEPIX[1]*page.rect.height,box[2]/PAGEPIX[0]*page.rect.width,box[3]/PAGEPIX[1]*page.rect.height)
    assert abs(rect.width-panel[0].rect.width)<.05 and abs(rect.height-panel[0].rect.height)<.05
    page.show_pdf_page(rect,panel,0,keep_proportion=False)
    prefix='Figure1E_original_canvas_position_proof' if number==1 else 'Figure3_COLI_SOX9_updated_A4'
    doc.save(ROOT/'figures'/f'{prefix}.pdf');doc.close()
    pix=panel[0].get_pixmap(matrix=fitz.Matrix((box[2]-box[0])/panel[0].rect.width,(box[3]-box[1])/panel[0].rect.height),alpha=False)
    image=Image.frombytes('RGB',(pix.width,pix.height),pix.samples).resize((box[2]-box[0],box[3]-box[1]),Image.Resampling.LANCZOS)
    proof=base.copy();proof.paste(image,(box[0],box[1]));mask=np.ones(PAGEPIX[::-1],bool);mask[box[1]:box[3],box[0]:box[2]]=False
    assert np.array_equal(original[mask],np.asarray(proof)[mask])
    proof.save(ROOT/'figures'/f'{prefix}.png',dpi=(600,600));proof.save(ROOT/'figures'/f'{prefix}.tiff',dpi=(600,600),compression='tiff_lzw')
    preview=proof.copy();preview.thumbnail((1800,2500));preview.save(ROOT/'audit'/f'{prefix}_preview.png')
    return {'figure':number,'slot_fraction':SLOTS[number],'slot_pixels':list(box),'slot_mm':[rect.width/72*25.4,rect.height/72*25.4],'outside_slot_pixel_changes':0,'base_canvas':str(SOURCE/'Figure 1.tif') if number==1 else str(CURRENT3),'scope':'Figure1E placement proof only; other Figure1 panels inherited from original TIFF, not current Figure1b-c' if number==1 else 'Current Figure3; only g replaced'}


def report(composites):
    counts=pd.concat([pd.read_csv(ROOT/'tables/Figure1E_main_sections_counts.csv'),pd.read_csv(ROOT/'tables/Figure3G_main_sections_counts.csv')],ignore_index=True)
    rows=[]
    for _,r in counts.iterrows():rows.append({'物种':r.species,'时间/阶段':r['section'] if r.species=='Human' else r.stage+' '+r.wound_time,'Fib空间观测':int(r.source_Fib_observations),'双阳空间观测':int(r.double_positive_Fib),'ROI内双阳':int(r.target_in_registered_zoom)})
    (ROOT/'tables/Current_panel_dual_positive_summary.csv').write_text(pd.DataFrame(rows).to_csv(index=False),encoding='utf-8-sig')
    # Copy the already-approved b-c unchanged for an explicit current-data reference, without shrinking 8 pt text.
    import shutil
    shutil.copy2(BC,ROOT/'figures/Figure1bc_current_unchanged_reference.png')
    def links(stem):return ' · '.join(f'<a href="../figures/{stem}.{ext}">{ext.upper()}</a>' for ext in ['pdf','svg','png','tiff'])
    panels=''
    for label,stem in [('Figure 1e（原占位大小）','Figure1E_COLI_SOX9_unified_slot'),('Figure 3g（原占位大小）','Figure3G_COLI_SOX9_unified_slot')]:
        panels+=f'<section><h2>{label}</h2><img src="../figures/{stem}.png"><p>{links(stem)}</p></section>'
    table=pd.DataFrame(rows).to_html(index=False,escape=True)
    body=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Figure 1e / 3g — COL I + SOX9</title><style>body{{margin:24px;background:#f3f6f8;color:#203848;font:15px Arial,"Microsoft YaHei",sans-serif}}section{{background:white;padding:24px;margin:18px auto;max-width:1240px}}h1{{font-size:26px}}h2{{font-size:22px}}img{{display:block;max-width:100%;height:auto}}.warn{{border-left:5px solid #bc9418}}table{{border-collapse:collapse}}td,th{{padding:8px 16px;border-bottom:1px solid #ddd}}.swatch{{display:inline-block;width:14px;height:14px;vertical-align:middle;margin-right:4px}}a{{color:#066389}}</style><h1>Figure 1e / Figure 3g：COL I＋SOX9 双marker统一版</h1>
    <section class="warn"><p><span class="swatch" style="background:{RED}"></span>SOX9 红色；<span class="swatch" style="background:{GREEN}"></span>COL I 绿色；<span class="swatch" style="background:{YELLOW}"></span>双阳黄色；其他空间点浅灰。原始免疫荧光照片未改色、未重绘。</p><p>空间图检测的是 RNA：人 <i>SOX9</i> &gt;0 且 <i>COL1A1</i> &gt;0；大鼠 <i>Sox9</i> &gt;0 且 <i>Col1a1</i> &gt;0，均限于原注释 Fib_ 谱系。COL1A1 是与 COL I 染色对应的转录本，不能将其等同于实测蛋白荧光强度。</p><p>原坐标、既有ROI、0.5 µm/坐标单位和比例尺沿用且不为双阳富集重新选择。Arial 8–12 pt；本次遵循原主图实际占位，独立面板PDF/SVG不额外缩小。虚线方框为放大ROI，不是瘢痕边界。</p><p>Figure1e主图仍显示2 mph与6 mph；6 mph两部分按旧主图几何合并，未当作两个独立样本。全部四份人输入（2 mph、6 mph两部分、9 mph）另有核对图。S/SL指已有锚点划定的scar/ scarless侧，不等同于独立健康对照。患者/动物独立编号不足，因此不做bin级伪重复检验：P=NE，FDR=NE。当前比例图为描述性结果。</p></section>
    {panels}<section><h2>当前两张面板的真实计数</h2>{table}<p><a href="../tables/Current_panel_dual_positive_summary.csv">计数CSV</a> · <a href="../tables/Human_four_source_parts_marker_audit.csv">人四份输入marker核对</a> · <a href="../tables/Figure3G_main_sections_counts.csv">大鼠切片metadata</a> · <a href="../tables/Human_marker_by_physical_section_region.csv">人切片区域比例（6 mph已合并）</a></p></section>
    <section><h2>双色marker核对：全部真实输入</h2><p>以下是转录本检测位置，不是生成的免疫染色或人工增加的共表达。</p><img src="../figures/Human_dual_marker_channel_audit_A4.png"><p>{links('Human_dual_marker_channel_audit_A4')}</p><img src="../figures/Rat_dual_marker_channel_audit_A4.png"><p>{links('Rat_dual_marker_channel_audit_A4')}</p></section>
    <section><h2>Figure 3整页A4：仅g替换，其他像素不变</h2><img src="../audit/Figure3_COLI_SOX9_updated_A4_preview.png"><p><a href="../figures/Figure3_COLI_SOX9_updated_A4.pdf">A4 PDF</a> · <a href="../figures/Figure3_COLI_SOX9_updated_A4.png">PNG 600 dpi</a> · <a href="../figures/Figure3_COLI_SOX9_updated_A4.tiff">TIFF</a></p></section>
    <section class="warn"><h2>Figure 1e原画布位置核对（不是新整张Figure 1）</h2><p>下面仅用于确认e的大小和位置，其余实验图片保持原TIFF；其中原TIFF的b–c仍是旧结果。已确认的新版b–c另列下方且文件完全不变，未缩小后强塞进旧槽位。请在最终主图源文件中使用独立e面板与新版b–c，勿将此位置核对页作为新统计整图投稿。</p><img src="../audit/Figure1E_original_canvas_position_proof_preview.png"><p><a href="../figures/Figure1E_original_canvas_position_proof.pdf">位置核对PDF</a></p><h2>已确认新版Figure 1b–c（原文件不变）</h2><img src="../figures/Figure1bc_current_unchanged_reference.png"></section>
    <section><h2>可复核性</h2><p><a href="../audit/INPUT_AND_MARKER_AUDIT.json">输入与marker审计</a> · <a href="../audit/RENDER_AND_CANVAS_QA.json">画布与绘图核查</a> · <a href="../outputs/SOX9_Figure1E_3G_COLI_SOX9_unified_v01.zip">全部结果ZIP</a></p><p>Table S1 单细胞清单未修改，本次空间数据metadata单独提供。原始数据及上一版结果未覆盖。</p></section></html>'''
    (ROOT/'06_html/OPEN_ME_FIRST_Figure1E_3G_COLI_SOX9_unified.html').write_text(body,encoding='utf8')
    qa={'status':'RENDERED_PENDING_FINAL_FILE_AND_VISUAL_QA','markers_only':['SOX9','COL1A1'],'cutoff':'both >0 inside source Fib_','palette':{'SOX9':RED,'COL_I_transcript':GREEN,'double_positive':YELLOW,'other':GRAY},'user_yellow_overlap_choice':True,'RNA_not_IF_measurement':True,'raw_IF_unmodified':True,'reference_functions_unchanged':True,'reference_sha256':sha(REFERENCE),'marker_and_coordinate_artifacts':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in (ROOT/'data').glob('*.npz')],'figures':OUTPUTS,'scale_bars_and_masks':BARS,'composites':composites,'current_Figure1bc_copy_sha256':sha(ROOT/'figures/Figure1bc_current_unchanged_reference.png'),'source_Figure1bc_sha256':sha(BC),'TableS1_changed':False,'stats':'Descriptive proportions; no independent patient/animal tests. P=NE and FDR=NE; 6 mph source parts pooled.'}
    (ROOT/'audit/RENDER_AND_CANVAS_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')


def main():
    regions=human_regions();base1,box1=human_main(regions);base3,box3=rat_main();channel_audits()
    composites=[compose(1,base1,box1,'Figure1E_COLI_SOX9_unified_slot'),compose(3,base3,box3,'Figure3G_COLI_SOX9_unified_slot')]
    report(composites);print('PANELS AND CANVAS PROOFS COMPLETE',flush=True)


if __name__=='__main__':main()
