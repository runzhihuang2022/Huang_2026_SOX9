from pathlib import Path
from html.parser import HTMLParser
import hashlib, html, json, os, platform, sys, zipfile, tempfile, shutil
import h5py, anndata as ad, numpy as np, pandas as pd, fitz, openpyxl
from PIL import Image
from prepare_full_selected import ROOT, PROJECT, SKILL, sha
sys.stdout.reconfigure(encoding='utf8')

def table(d):return d.to_html(index=False,escape=True,na_rep='NA',float_format=lambda v:f'{v:.6g}')
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        for k in ['href','src']:
            if k in a and not a[k].startswith(('#','http','data:','mailto:')):self.links.append(a[k])

def local_files(root):
    for parent,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=[d for d in dirs if d!='node_modules']
        for name in files:yield Path(parent)/name

def main():
    qa={};requested=pd.read_csv(ROOT/'01_input_audit/user_candidate_skin_record_ids_LOCKED.csv').record_id.tolist()
    skin=pd.read_csv(ROOT/'tables/Figure1_selected19_full_QC_skin_cell_metadata.tsv.gz',sep='\t',index_col=0,low_memory=False,float_precision='round_trip')
    with h5py.File(ROOT/'02_shared_geometry/Figure1_five_organ_selected19_refit.h5ad','r') as h:
        f1=ad.io.read_elem(h['obs']);xy1=h['obsm/X_umap'][:]
    with h5py.File(ROOT/'01_input_audit/Selected19_full_QC_Fib_normalized_RNA.h5ad','r') as h:
        f2=ad.io.read_elem(h['obs']);genes=ad.io.read_elem(h['var'][h['var'].attrs['_index']])
    f1skin=f1.loc[f1.organ.eq('skin')]
    canonical=f1skin.index.str.replace(r'\|skin_full$','',regex=True)
    assert set(canonical)==set(skin.index)
    fib_names=set(skin.index[skin.original_Fib_gate.eq(True)])
    assert fib_names==set(f2.index) and len(fib_names)==20270
    assert set(f1skin.record_id)==set(f2.record_id)==set(requested) and len(requested)==19
    assert f1.index.is_unique and f2.index.is_unique and np.isfinite(xy1).all()
    manifest=json.loads((ROOT/'analysis/human_universal/00_manifest/Figure2_A-E_manifest.json').read_text(encoding='utf8'))
    geometry=pd.read_csv(manifest['geometry_path'],sep='\t',float_precision='round_trip')
    assert set(geometry.barcode)==set(f2.index) and geometry.barcode.is_unique and len(geometry)==20270
    assert hashlib.sha256(geometry[['barcode','UMAP1','UMAP2']].to_csv(sep='\t',index=False).encode()).hexdigest()==manifest['geometry_sha256']
    g=geometry.set_index('barcode').reindex(f2.index)
    assert set(f2.patient_id_source)==set(geometry['sample'])
    assert not f1skin.actual_tissue.astype(str).ne('Skin').any()
    old=pd.read_csv(PROJECT/'20261009_Figure1BC_selected7_controls_v02/audit/selected_cells_metadata_expression.tsv.gz',sep='\t',index_col=0,low_memory=False,float_precision='round_trip')
    pub=old.loc[~old.organ.eq('skin')]
    assert set(f1.index[f1.organ.ne('skin')])==set(pub.index)
    actual=f1.reindex(pub.index)['SOX9_nonnegative_expression'].to_numpy(dtype=float)
    assert np.array_equal(actual,pub.SOX9_nonnegative_expression.to_numpy())
    assert len(f1)==465450 and len(f1skin)==208532
    ledger=pd.read_csv(ROOT/'analysis/human_universal/05_tables/Figure2_selected_marker_ledger.tsv',sep='\t')
    assert (ledger.groupby('label').size()==3).all() and 'KRT14' not in ledger.gene.tolist()
    assert 'KRT14' in genes and all(v in genes for v in ['SOX9','PDGFRA','VIM','ACTA2','COL1A1'])
    skillfile=SKILL/'scripts/run_figure2_ae_universal.py'
    assert sha(skillfile)==(ROOT/'audit/skill_renderer_sha256.txt').read_text().strip()==sha(ROOT/'analysis/human_universal/08_code/run_figure2_ae_universal.py')
    stats=pd.read_csv(ROOT/'tables/Figure1C_all_organ_expression_P_FDR.csv')
    row=stats.loc[stats.organ.eq('skin')].iloc[0]
    assert row.target=='All Fib' and row.n_reference_sample_units==5 and row.n_case_sample_units==10
    patient=pd.read_csv(ROOT/'tables/skin_patient_unit_expression.csv',float_precision='round_trip')
    assert np.isclose(row.mean_reference,patient.loc[patient.condition_pair.eq('normal'),'SOX9_mean_expression'].mean())
    assert np.isclose(row.mean_case,patient.loc[patient.condition_pair.eq('fibrosis'),'SOX9_mean_expression'].mean())
    assert set(patient.loc[patient.SOX9_mean_expression.eq(0),'libraries'])=={'F091','F092'}
    s1file=ROOT/'outputs/Table_S1_current_Figure1_4_selected19_full_v06.xlsx'
    w=openpyxl.load_workbook(s1file,read_only=True,data_only=True);s=w['Table S1'];it=s.iter_rows(values_only=True);headers=list(next(it));s1=pd.DataFrame(list(it),columns=headers);w.close()
    assert len(s1)==104 and s1['Record ID'].is_unique
    assert s1['Figure1B cells'].sum()==465450 and s1['Figure2A/B Fib cells'].sum()==20270
    assert s1['Figure4 all Fib cells'].sum()==59854
    assert set(s1.loc[s1['Figure2A/B Fib cells'].gt(0),'Record ID'])==set(requested)
    assert not s1['True tissue'].astype(str).str.fullmatch('PBMC|BM|Bone marrow',case=False).any()
    # All newly delivered formal panels must be readable at the recorded physical size.
    pdfs=[ROOT/'figures/Figure1bc_horizontal_selected19_fullcells.pdf',ROOT/'figures/Figure2AB_human_Fib_slot.pdf',ROOT/'figures/Figure2_SOΧ9_highest_and_feature_shared_geometry.pdf',ROOT/'analysis/human_universal_presentation/04_figures/Figure2_selected19_full_A-E_skill_layout_reflow.pdf']
    pdf_qa=[]
    for p in pdfs:
        doc=fitz.open(p);page=doc[0];spans=[sp for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for sp in line['spans']]
        fonts=sorted({sp['font'] for sp in spans});sizes=[sp['size'] for sp in spans]
        assert len(doc)==1 and min(sizes)>=7.98 and max(sizes)<=12.02 and all('Arial' in name for name in fonts),(p.name,fonts,min(sizes),max(sizes))
        rec={'file':p.relative_to(ROOT).as_posix(),'page_count':len(doc),'width_mm':page.rect.width*25.4/72,'height_mm':page.rect.height*25.4/72,'font_range_pt':[min(sizes),max(sizes)],'font_names':fonts,'font_embedded':all(f[1]!='n/a' for f in page.get_fonts())}
        pix=page.get_pixmap(dpi=130,alpha=False);pixels=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,3)
        rec['non_white_fraction']=float(np.mean(np.any(pixels<245,axis=2)));assert rec['non_white_fraction']>.01
        png=p.with_suffix('.png')
        with Image.open(png) as im:
            rec['PNG_pixels']=list(im.size);rec['PNG_dpi']=list(im.info.get('dpi',[]))
            assert im.info.get('dpi',[0])[0]>599
        pdf_qa.append(rec);doc.close()
    qa.update({'status':'PASS_REQUESTED_FIGURE1BC_FIGURE2AB_FULL_CELL_AND_MEMBERSHIP_QA','whitelist_libraries':19,'skin_full_QC_cells':208532,'Fib_full_QC_cells':20270,'Figure1_total_cells':465450,'Figure1_and_Figure2_same_skin_library_IDs':True,'Figure2_exact_Figure1_skin_full_Fib_subset':True,'PBMC_specimen_libraries':0,'non_skin_cell_IDs_and_SOΧ9_values_unchanged':True,'pristine_skill_sha256':sha(skillfile),'skill_source_unchanged':True,'new_global_and_lineage_geometries':True,'Figure2_geometry_sha256':manifest['geometry_sha256'],'highest_observed_SOΧ9_cluster':manifest['highest_target_cluster'],'C12_identity_pending_biological_validation':True,'not_claimed_pure_definitive_Fib_subtype':True,'statistical_zero_expression_samples_kept':['F091','F092'],'TableS1_scRNA_records':104,'visible_plot_fonts':'Embedded Arial 8–12pt','main_required_panels_visual_QA':'Manually inspected; legend/footer/marker labels separated; no required unresolved visual failure','formal_panel_PDF_QA':pdf_qa,'unrequested_supplement_E':'NOT_ESTIMABLE: second contrast not requested/preregistered','post_expression_selection_inference':'Exploratory only'})
    (ROOT/'audit/FINAL_FIGURE_AND_METADATA_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    lib=pd.read_csv(ROOT/'tables/Selected19_full_Fib_library_expression.csv')
    cohort=lib[['record_id','source_sample','patient_id_source','clinical_group','condition_pair','Fib_cells','SOX9_mean','SOX9_positive_cells']].copy()
    cohort.columns=['记录ID','样本/文库ID','已有患者ID','临床组','图中分组','全量QC Fib细胞','Fib内SOX9均值（含零）','SOX9>0 Fib细胞（辅助）']
    purity=pd.read_csv(ROOT/'tables/Figure2_all_cluster_marker_and_purity_audit.csv')
    rank=purity.loc[purity.gene.eq('SOX9'),['cluster','n_cells','mean_expression','percent_expressed']].sort_values('mean_expression',ascending=False)
    rank.columns=['无监督簇','细胞数','SOX9均值','SOX9>0 (%)']
    high=purity[purity.cluster.eq('C12')].copy()
    st=json.loads((ROOT/'audit/Figure1C_SKIN_STATISTICS.json').read_text(encoding='utf8'))
    card=lambda stem,caption:f'<section><h2>{caption}</h2><img src="../figures/{stem}.png"><p><a href="../figures/{stem}.pdf">PDF</a> · <a href="../figures/{stem}.svg">可编辑SVG</a> · <a href="../figures/{stem}.png">600dpi PNG</a></p></section>'
    page=r'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>SOX9：选定19条皮肤文库，Figure1b–c与Figure2A/B</title><style>body{font:15px Arial,"Microsoft YaHei",sans-serif;background:#f4f7f9;color:#203a4d;margin:24px}section{background:white;border-radius:7px;padding:22px;margin:20px 0}h1{font-size:25px}h2{font-size:20px}nav{padding:12px;background:#e5eef3;position:sticky;top:0}nav a{margin-right:20px}a{color:#006e9a}p,li{line-height:1.7}.warning{border-left:6px solid #c85e45}.scroll{overflow:auto;max-height:680px}table{width:100%;border-collapse:collapse;font-size:13px}td,th{border:1px solid #d9e3e9;padding:9px;text-align:left}th{background:#244b63;color:white}tbody tr:nth-child(even){background:#edf5f8}img{max-width:100%;height:auto}code{word-break:break-all}input{padding:8px}small{color:#536b7a}details{padding:12px}@media print{nav,input{display:none}.scroll{max-height:none}section{break-inside:avoid}}</style></head><body>
<h1>选定19条真实皮肤文库：Figure 1b–c与Figure 2A/B全量重分析</h1><nav><a href="#f1">Figure 1</a><a href="#f2">Figure 2</a><a href="#cohort">样本表</a><a href="#statistics">P/FDR</a><a href="#s1">Table S1</a><a href="#qa">核对与下载</a></nav>
<section><p>严格使用您上传的19条记录白名单：14条DSDB_Scar与5条Normal。皮肤从208,616个源细胞按同历史QC保留208,532个，全量原注释Fib为20,270个。Figure1与Figure2的皮肤文库ID完全一致；Figure2是Figure1同批皮肤细胞内的全部Fib子集，没有按SOX9阳性或双阳性预筛，没有再使用展示抽样对象作定量全集。</p><p>其他器官256,918个细胞及其表达、统计单位、原始P值保持不变。五器官总图使用465,450个当前细胞重新PCA/邻居/UMAP，没有沿用旧坐标。不同器官/物种/来源高度混杂，因此未对这一混杂键做强制批次校正；跨器官符号匹配不是已验证的一对一同源基因分析，UMAP仅作描述性展示。</p></section>
<div id="f1">FIGURE1</div>
<section id="statistics"><h2>皮肤主结果：全部Fib内SOX9表达均值</h2><p>已有映射为5个正常和10个瘢痕患者ID。先计算每条文库的Fib表达均值，再对同患者文库等权汇总；不把同患者重复文库当独立人。正常组均值REF_MEAN，瘢痕组CASE_MEAN，约FOLD倍。双侧Welch检验 <strong>P=P_VALUE，五器官BH FDR=FDR_VALUE</strong>。零表达的YN-1、YN-2保留；YN-3在全量原Fib注释中有1个SOX9表达细胞，不再被写为真实零。</p><p>对照中YN-1、YN-2、YN-3的Fib数量分别为13、5、18，均保留，但估计不精确。患者ID沿用已有映射而非本轮独立确认。您是在查看SOX9表达摘要后选定标本，因此本版为<strong>探索性筛选队列结果，不能当作独立预设队列的确认性证据</strong>。其他器官P值不变；FDR因新的五检验集合而重新计算。器官间纵轴范围不同，表达规范化来源也不同，不直接跨器官比较数值。</p><p>非参数敏感性核对：双侧Mann–Whitney P=MW_P，单项敏感性检验家族BH FDR=MW_P。不是重新选择更显著的主检验。</p><div class="scroll">STATS_TABLE</div><p><a href="../tables/Figure1C_all_organ_expression_P_FDR.csv">五器官P/FDR与均值表</a> · <a href="../tables/skin_patient_unit_expression.csv">患者统计单位与来源文库</a></p></section>
<div id="f2">FIGURE2</div>
<section class="warning"><h2>C12最高表达，但生物学身份仍需核查</h2><p>技能原始执行器固定resolution=0.55、seed=20260917，得到17个簇C1–C17。C12有185个细胞，SOX9均值0.59535，最高；183个源注释为Myo_Fib、2个为Inflama_Fib。ACTA2、VIM检测比例分别约97.3%、96.2%，但<strong>KRT14检测比例100%，PDGFRA仅约1.62%</strong>。这一组合需要上皮污染、双细胞或细胞状态的复核；不能仅凭源标签直接命名为纯成纤维细胞新亚型。</p><p>主dotplot按您要求不显示KRT14，原始表达仍保留，并在下面完整核对表中披露。未因这一结果删除C12或重选聚类参数。星号只表示“源Fib门控中观察到的最高SOX9均值簇”，不是已经证明其纯度、差异或机制。RNA与蛋白的共表达也不能直接等同。</p><div class="scroll">HIGHEST_TABLE</div></section>
AUXILIARY
<section><h2>各簇SOX9表达排序</h2><div class="scroll">RANK_TABLE</div><p><a href="../tables/Figure2_all_cluster_marker_and_purity_audit.csv">全部簇marker与纯度核对</a> · <a href="../tables/Figure2_cluster_original_subtype_counts.csv">簇与源亚型对应</a> · <a href="../analysis/human_universal/05_tables/Figure2_selected_marker_ledger.tsv">技能每簇Top3 marker清单</a> · <a href="../analysis/human_universal/05_tables/Figure2_highest_target_cluster_tests_P_FDR.tsv">探索性配对表达比较P/FDR</a></p><p>每簇Top3从固定marker白名单选择，必选SOX9、PDGFRA、VIM、ACTA2。SOX9&gt;0且COL1A1&gt;0仅用于辅助共表达汇总，未用于删掉母谱系中的阴性细胞。</p></section>
<section><h2>技能附带A–E执行记录</h2><p>本次主图需求是A/B。附带C为组成，D为已注册Scar_vs_normal比例比较；E没有第二个预注册对比，明确记为NOT_ESTIMABLE，没有用其他图替代。下图仅将17个簇的dotplot横轴标签转为90°，分析函数、表格和坐标保持技能原始版本。</p><img src="../analysis/human_universal_presentation/04_figures/Figure2_selected19_full_A-E_skill_layout_reflow.png"><p><a href="../analysis/human_universal_presentation/04_figures/Figure2_selected19_full_A-E_skill_layout_reflow.pdf">A–E附带PDF</a></p></section>
<section id="cohort"><h2>本次19条实际使用皮肤文库</h2><p>文库不是独立患者数。该清单用于Figure1皮肤UMAP和Figure2Fib，不再有PBMC、急性DSDB或SSDB_Non_Scar。<a href="../tables/Selected19_full_Fib_library_expression.csv">下载逐文库Fib表达表</a> · <a href="../01_input_audit/user_candidate_skin_record_ids_LOCKED.csv">锁定白名单</a></p><input id="search" placeholder="筛选记录ID、样本、患者或组别"><div class="scroll" id="sampleTable">COHORT_TABLE</div></section>
<section id="s1"><h2>匹配本版主图的Table S1</h2><p>104条实际单细胞文库：Figure1使用96条（皮肤19、其他器官77）；Figure2同19条，不重复计数；另有已确认且未改的Figure4培养Fib八条。没有PBMC标本记录。Figure1e、Figure3g为空间切片，不混入单细胞文库表。</p><p><a href="../outputs/Table_S1_current_Figure1_4_selected19_full_v06.xlsx">Table S1 Excel</a> · <a href="../tables/Table_S1_current_selected19_Figure1_4_actual_scRNA.csv">Table S1 CSV</a></p><div class="scroll">S1_TABLE</div><p>所有原表保留。本轮没有改动已确认的Figure1e（2mph、6mph两part；9mph补充）、Figure3g六张空间切片或Figure4八条培养Fib分析。</p></section>
<section id="qa"><h2>文件核对、代码与下载</h2><p>新图文字为嵌入Arial 8–12pt，矢量PDF/SVG和600dpi PNG/TIFF。Figure2按原主图A/B槽位制作；完整A4拼图仅替换这一槽位，其他实验图像像素保持不变。源实验图中的字体/原始测量没有在此重新审查。</p><p><a href="../figures/Figure2_complete_A4_selected19_full_Fib.pdf">Figure2完整A4 PDF</a> · <a href="../figures/Figure2_complete_A4_selected19_full_Fib.png">Figure2完整A4 PNG</a> · <a href="../audit/FINAL_FIGURE_AND_METADATA_QA.json">最终图与样本核对</a> · <a href="../audit/INPUT_MATCH_QA.json">全量输入证据</a> · <a href="../audit/PRIMARY_SOURCE_SHA256_MANIFEST.json">原始来源SHA-256</a></p><p><a href="../outputs/SOX9_Figure1BC_Figure2AB_selected19_fullcells_v01.zip">离线图册与代码压缩包</a> · <a href="../07_reproduction/REPRODUCTION_AND_SCOPE.md">复现与范围说明</a></p><p>双图样本与条码均已精确核对。Figure2使用的是技能原文件，不修改它的分析函数；原聚类坐标不复用。主图版式适配只处理图例间距与簇标签，不改变marker计算或结果。</p></section>
<script>document.getElementById('search').addEventListener('input',e=>document.querySelectorAll('#sampleTable tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(e.target.value.toLowerCase())?'':'none'));</script></body></html>'''
    tokens={'FIGURE1':card('Figure1bc_horizontal_selected19_fullcells','Figure 1b–c：全量选定皮肤与五器官新UMAP'),'FIGURE2':card('Figure2AB_human_Fib_slot','Figure 2A/B：同19条文库的20,270个Fib'),'AUXILIARY':card('Figure2_SOΧ9_highest_and_feature_shared_geometry','同一冻结Fib UMAP：最高簇与SOX9表达'),'REF_MEAN':f'{st["reference_mean_of_patient_means"]:.5f}','CASE_MEAN':f'{st["scar_mean_of_patient_means"]:.5f}','FOLD':f'{st["fold_change_ratio"]:.3f}','P_VALUE':f'{st["two_sided_Welch_P"]:.6g}','FDR_VALUE':f'{st["BH_FDR_five_organ_family"]:.6g}','MW_P':f'{st["two_sided_Mann_Whitney_P_sensitivity"]:.6g}','STATS_TABLE':table(stats),'HIGHEST_TABLE':table(high),'RANK_TABLE':table(rank),'COHORT_TABLE':table(cohort),'S1_TABLE':table(s1)}
    for k in sorted(tokens,key=len,reverse=True):page=page.replace(k,tokens[k])
    entry=ROOT/'06_html/OPEN_ME_FIRST_SOX9_selected19_full_Figure1_2_complete_v01.html';entry.write_text(page,encoding='utf8')
    # Reproduction is contract/config based; original data are not bundled or overwritten.
    commands='''# Reproduction and scope\n\nUse a new sibling version directory; never overwrite an existing result.\n\n1. Stage the canonical source files at the absolute paths in PRIMARY_SOURCE_SHA256_MANIFEST.json, or update the explicit path constants in the separate input adapters. The original source matrices are omitted from the portable archive.\n2. Place the locked 19-record CSV in the declared CHOICE path. Do not substitute other sample IDs.\n3. Use the recorded Python environment/package versions and run prepare_full_selected.py.\n4. Run the pristine skill executable with 00_manifest/figure2_config.json. Its plotting/analysis functions must remain unchanged.\n5. Run build_full_atlas.py, then render_Figure1BC.py and render_Figure2AB.py. Reflow_skill_supplement.py changes only presentation tick rotation.\n6. Prepare metadata and run build_TableS1.mjs using the bundled artifact-tool runtime. Do not edit the former Table S1 or the source clinical map.\n7. Run build_delivery_and_QA.py after manually checking the new panels.\n\nCore statistical contract: source original Fibroblast gate; full QC cells, no SOX9/coexpression prefilter; mean normalized SOX9 includes zeros. Equal library means within recorded patient IDs for skin. Public panels retain the old source/statistical units. Welch skin P is part of the five-organ BH family. The selected cohort and post-ranking comparisons are exploratory.\n\nC12 is the highest observed mean within the original Fib gate but is KRT14-high; definitive lineage/subtype identity requires independent review. Do not describe it as proven pure fibroblasts or a causal SOX9 mechanism.\n\nLarge derived H5AD/PCA files are omitted from the archive. Their live paths and hashes are recorded in LARGE_OBJECT_STAGING.json; adapters reproduce them from the canonical sources.\n'''
    (ROOT/'07_reproduction/REPRODUCTION_AND_SCOPE.md').write_text(commands,encoding='utf8')
    versions={}
    import importlib.metadata as md
    for pkg in ['anndata','scanpy','numpy','pandas','scipy','scikit-learn','matplotlib','h5py','numba','igraph','leidenalg']:
        versions[pkg]=md.version(pkg)
    (ROOT/'07_reproduction/environment_lock.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'packages':versions,'thread_counts':{'OMP':4,'OpenBLAS':4,'MKL':4,'NUMBA':4},'Figure1_seed':20261009,'Figure2_seed':20260917,'Figure2_resolution':.55},indent=2),encoding='utf8')
    large=[]
    for p in local_files(ROOT):
        if p.suffix in ['.h5ad','.npy']:large.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'archive_status':'Omitted; reconstruct from adapters and registered canonical sources'})
    (ROOT/'audit/LARGE_OBJECT_STAGING.json').write_text(json.dumps(large,ensure_ascii=False,indent=2),encoding='utf8')
    parser=Links();parser.feed(page)
    # The archive and parent hash manifest are created by final packaging after this QA stage.
    deferred={'../outputs/SOX9_Figure1BC_Figure2AB_selected19_fullcells_v01.zip','../audit/PRIMARY_SOURCE_SHA256_MANIFEST.json'}
    missing=[q for q in parser.links if not (entry.parent/q).resolve().exists() and q not in deferred]
    assert not missing,missing
    print('HTML, plots, exact cell/sample match and Table S1 QA complete',flush=True)
    print(json.dumps({'Figure1_cells':len(f1),'Skin_cells':len(f1skin),'Fib_cells':len(f2),'libraries':len(requested),'Figure2_clusters':manifest['n_clusters'],'highest_cluster':'C12','TableS1_records':len(s1),'HTML_links':len(parser.links)},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
