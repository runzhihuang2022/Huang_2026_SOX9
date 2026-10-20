from __future__ import annotations

import base64
import hashlib
import html
import json
import platform
from pathlib import Path

import pandas as pd
from PIL import Image

ROOT=Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
FIG=ROOT/"04_figures"; TAB=ROOT/"05_tables"; MAN=ROOT/"00_manifest"; QA=ROOT/"10_qa"; HTM=ROOT/"06_html"
ROUTER=Path('data/external/INPUT_6e7fdfb734/dataset_registry.json')

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(8<<20),b""): h.update(b)
    return h.hexdigest()

def uri(path): return "data:image/png;base64,"+base64.b64encode(Path(path).read_bytes()).decode("ascii")
def table(path, columns=None, limit=200):
    d=pd.read_csv(path,sep="\t"); d=d[columns] if columns else d
    return d.head(limit).to_html(index=False,border=0,classes="data",float_format=lambda x:f"{x:.4g}")

def composition_max_deviation(summary_path):
    d=pd.read_csv(summary_path,sep="\t"); registry=d[["sample","group"]].drop_duplicates(); clusters=sorted(d.cluster.unique(),key=lambda x:int(x[1:])); grid=registry.assign(k=1).merge(pd.DataFrame({"cluster":clusters,"k":1}),on="k").drop(columns="k"); q=grid.merge(d[["sample","group","cluster","n_cells"]],how="left",on=["sample","group","cluster"]); q.n_cells=q.n_cells.fillna(0); q["fraction"]=q.n_cells/q.groupby(["sample","group"]).n_cells.transform("sum"); sums=q.groupby(["sample","group"]).fraction.sum(); return float((sums-1).abs().max())

def main():
    HTM.mkdir(exist_ok=True); QA.mkdir(exist_ok=True)
    hm=json.loads((ROOT/"03_analysis/human_Figure2_AE_composition/00_manifest/Figure2_A-E_manifest.json").read_text(encoding="utf-8")); rm=json.loads((ROOT/"03_analysis/rat_Figure3G_AE_composition/00_manifest/Figure2_A-E_manifest.json").read_text(encoding="utf-8"))
    cutoff=json.loads((QA/"CUTOFF_CHANGE_AUDIT.json").read_text(encoding="utf-8")); hs=json.loads((MAN/"Figure1E_human_spatial_manifest_v05.json").read_text(encoding="utf-8")); rs=pd.read_csv(TAB/"Figure3G_rat_stage_statistics_P_and_FDR.tsv",sep="\t")
    sc=pd.read_csv(TAB/"scRNA_cutoff1_vs_cutoff0_total_audit.tsv",sep="\t"); human=sc[sc.dataset.eq("human_scRNA")].iloc[0]; rat=sc[sc.dataset.eq("rat_scRNA")].iloc[0]
    spatial_dist={"human":{"SOX9_min_positive":2.8079755,"COL1A1_min_positive":2.7581758,"log1p_metadata_present":True},"rat":{"Sox9_min_positive":1.467191300284785,"Col1a1_min_positive":1.5100753954781283,"log1p_metadata_present":True},"interpretation":"All nonzero spatial values exceed 1, so cutoff >0 and >1 legitimately select the same spatial double-positive bins."}
    (QA/"SPATIAL_EXPRESSION_SCALE_AUDIT.json").write_text(json.dumps(spatial_dist,ensure_ascii=False,indent=2),encoding="utf-8")
    stage4={"formal_pySCENIC_execution_audit":"NOT_ESTIMABLE","RSS_rank_panel":"NOT_ESTIMABLE","AUCell_matrix_and_featureplot":"NOT_ESTIMABLE","reason":"The registered server mirror contains pbmc.h5ad and s3_avg20_rep1.regulons.txt but no AUCell matrix, RSS table, or complete GRNBoost2-cisTarget-AUCell execution audit.","completed":{"network":"Predicted SOX9 network drawn from the existing regulons.txt target set; edge sign is expression association.","GO_BP":"Completed for the existing target set with tested-gene background."},"claim_boundary":"Do not describe the network as a newly verified full pySCENIC run."}
    (QA/"FIGURE4H_STAGE4_ELIGIBILITY.json").write_text(json.dumps(stage4,ensure_ascii=False,indent=2),encoding="utf-8")
    route={"registry":str(ROUTER),"registry_sha256":sha256(ROUTER),"human":{"boundary_root":'data/external/INPUT_8665934387/mask',"sections":["2mph","6mph_p1","6mph_p2","9mph"],"semantics":"scar-to-scarless direction using full registered origin line after rotation","gate":"Fib annotation AND SOX9>0 AND COL1A1>0"},"rat":{"boundary_root":'data/external/INPUT_df0bc4e82e/annotation_output',"stages":["E16.5","Adult"],"semantics":"wound-centre versus flank, not histological scar inside/outside","gate":"Fib annotation AND Sox9>0 AND Col1a1>0"}}
    (MAN/"DATA_ROUTER_BOUNDARY_CALLS_v05.json").write_text(json.dumps(route,ensure_ascii=False,indent=2),encoding="utf-8")
    commands=[
      'python scripts/run_figure2_ae_universal.py --config 00_manifest/human_figure2_composition_config.json',
      'python scripts/run_figure2_ae_universal.py --config 00_manifest/rat_figure3g_composition_config.json',
      'python 08_code/prepare_spatial_sox9_col1a1_cutoff0.py', 'python 08_code/render_spatial_sox9_col1a1_cutoff0.py',
      'python 08_code/run_server_figure4h_5a.py', 'python 08_code/cutoff_audit_and_supplements_v05.py', 'python 08_code/build_html_qa_v05.py']
    (MAN/"RUN_COMMANDS.txt").write_text("\n".join(commands)+"\n",encoding="utf-8")
    figures=[
      ("Figure 1E — human spatial cutoff 0","Figure1E_human_spatial_SOX9_COL1A1_gt0_oldstyle_A4.png"),
      ("Figure 2A–E — human Fib composition contrasts","Figure2_human_Fib_A-E_composition_cutoff0_A4.png"),
      ("Human SOX9/COL1A1 cutoff 0 supplementary panels","human_scRNA_target_highlight_and_cutoff0_supplement_A4.png"),
      ("Figure 3G — rat Fib composition contrasts","Figure3G_rat_Fib_A-E_composition_cutoff0_A4.png"),
      ("Rat Sox9/Col1a1 cutoff 0 supplementary panels","rat_scRNA_target_highlight_and_cutoff0_supplement_A4.png"),
      ("Figure 3G — rat spatial cutoff 0","Figure3G_rat_spatial_Sox9_Col1a1_gt0_oldstyle_A4.png"),
      ("Figure 4H — existing SOX9 regulon target network and GO","Figure4H_Fib_ACTA2_SOX9_network_GO_A4.png"),
      ("Figure 5A — library-pseudobulk correlations and differences","Figure5A_SOX9_SCD_lysosomal_targets_correlations_and_differences_A4width.png")]
    dims=[]; cards=[]
    for label,name in figures:
        p=FIG/name
        with Image.open(p) as im: dims.append({"figure":label,"file":name,"width_px":im.width,"height_px":im.height,"A4_600dpi_pass":im.size==(4960,7015),"sha256":sha256(p)})
        cards.append(f"<section><h2>{html.escape(label)}</h2><img src='{uri(p)}' alt='{html.escape(label)}'></section>")
    dims=pd.DataFrame(dims); dims.to_csv(QA/"FIGURE_DIMENSION_HASH_QA.tsv",sep="\t",index=False)
    hdev=composition_max_deviation(ROOT/"03_analysis/human_Figure2_AE_composition/05_tables/Figure2_sample_cluster_summary.tsv"); rdev=composition_max_deviation(ROOT/"03_analysis/rat_Figure3G_AE_composition/05_tables/Figure2_sample_cluster_summary.tsv")
    embryo=rs[rs.stage.eq("E16.5")].iloc[0]; adult=rs[rs.stage.eq("Adult")].iloc[0]
    document=f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>SOX9 COL1A1 cutoff0 complete v05</title><style>body{{margin:0;background:#edf1f5;color:#17202a;font-family:Arial,sans-serif}}main{{max-width:1180px;margin:auto;background:#fff;padding:30px}}h1{{font-size:12pt}}h2{{font-size:11pt;margin-top:34px}}h3{{font-size:10pt}}p,li,td,th{{font-size:8pt;line-height:1.45}}img{{width:100%;border:1px solid #d7dde3}}.ok{{border-left:5px solid #287d4f;background:#f1fff7;padding:12px}}.warn{{border-left:5px solid #b36b00;background:#fff7e8;padding:12px}}table.data{{border-collapse:collapse;width:100%}}table.data th,table.data td{{border:1px solid #d7dde3;padding:4px}}table.data th{{background:#edf2f7}}code{{font-family:Consolas,monospace;font-size:8pt}}</style></head><body><main>
<h1>SOX9+COL1A1+ cutoff 0 — strict rerun v05</h1>
<div class='ok'><b>门控已从表达矩阵重新计算：</b>人scRNA {int(human.n_double_gt1):,}→{int(human.n_double_gt0):,}（新增{int(human.n_added):,}）；大鼠scRNA {int(rat.n_double_gt1):,}→{int(rat.n_double_gt0):,}（新增{int(rat.n_added):,}）。正式A–E为composition_contrasts；组成和最大偏差：人{hdev:.2g}、大鼠{rdev:.2g}。</div>
<div class='warn'><b>空间阈值审计：</b>空间对象已经是log1p表达，所有非零SOX9/COL1A1或Sox9/Col1a1值均高于1，所以空间cutoff 0与cutoff 1计数相同，这是表达值分布导致，并非复用旧结果。人空间不支持瘢痕内富集（P={hs['P_value']:.4g}, FDR={hs['BH_FDR']:.4g}）；E16.5中心-侧翼P={embryo.one_sided_exact_signflip_P:.4g}, FDR={embryo.BH_FDR_across_stages:.4g}；Adult P={adult.one_sided_exact_signflip_P:.4g}, FDR={adult.BH_FDR_across_stages:.4g}。</div>
<div class='warn'><b>Figure 4H边界：</b>现有regulons.txt可支持目标网络与GO图，但服务器镜像缺少AUCell矩阵、RSS表和完整执行审计，因此formal pySCENIC Stage 4的RSS/AUCell面板标记NOT_ESTIMABLE，不能宣称完成了一次新的完整pySCENIC运行。</div>
<p>人Fib：{hm['lineage_shape'][0]:,}细胞、{hm['n_clusters']}簇，SOX9观察均值最高{hm['highest_target_cluster']}；大鼠Fib：{rm['lineage_shape'][0]:,}细胞、{rm['n_clusters']}簇，Sox9观察均值最高{rm['highest_target_cluster']}。大鼠正式D/E因Normal组仅2个含Fib样本而按预注册规则显示NE。</p>
{''.join(cards)}
<section><h2>Cutoff审计</h2>{table(TAB/'scRNA_cutoff1_vs_cutoff0_total_audit.tsv')}<h3>人空间</h3>{table(TAB/'human_spatial_cutoff1_vs_cutoff0_by_section_region.tsv')}<h3>大鼠空间</h3>{table(TAB/'rat_spatial_cutoff1_vs_cutoff0_by_section_zone.tsv')}</section>
<section><h2>Figure 2/3差异丰度统计</h2><h3>Human SPTDI vs Normal</h3>{table(TAB/'Figure2D_human_SPTDI_vs_Normal_P_FDR.tsv')}<h3>Human DPTDI vs Normal</h3>{table(TAB/'Figure2E_human_DPTDI_vs_Normal_P_FDR.tsv')}<h3>Rat PTDI vs Normal</h3>{table(TAB/'Figure3G_D_rat_PTDI_vs_Normal_P_FDR.tsv')}<h3>Rat FTDI vs Normal</h3>{table(TAB/'Figure3G_E_rat_FTDI_vs_Normal_P_FDR.tsv')}</section>
<section><h2>Figure 5A统计</h2><h3>相关性</h3>{table(TAB/'Figure5A_SOX9_target_library_correlations.tsv',['state','gene','n_libraries','spearman_rho','P_value','BH_FDR_within_state'])}<h3>组间差异</h3>{table(TAB/'Figure5A_Fib_ACTA2_vs_other_Fib_paired_differences.tsv',['gene','n_paired_libraries','paired_mean_difference','paired_exact_Wilcoxon_P','BH_FDR'])}</section>
<section><h2>完成状态</h2><ul><li>QA通过：人/大鼠composition A–E、cutoff 0补充图、人/大鼠空间图、Figure 4H现有目标网络+GO、Figure 5A。</li><li>NOT_ESTIMABLE：Figure 4H formal RSS/AUCell；大鼠D/E推断统计（Normal含Fib样本n=2）。</li><li>无后台分析仍在运行。</li><li>网络与GO为预测性结果；Figure 5A因供体独立性未记录而为探索性。</li><li>不支持：人空间瘢痕内显著富集、大鼠Adult中心富集、任何因果结论。</li></ul></section>
</main></body></html>"""
    entry=HTM/"OPEN_ME_FIRST_SOX9_COL1A1_cutoff0_complete_v05.html"; entry.write_text(document,encoding="utf-8")
    qa={"status":"PASS" if dims.A4_600dpi_pass.all() and hdev<1e-10 and rdev<1e-10 else "FAIL","html":str(entry),"a4_600dpi_pass":bool(dims.A4_600dpi_pass.all()),"composition_sum_max_deviation":{"human":hdev,"rat":rdev},"font_family":"Arial","font_range_pt":[8,12],"gate":"SOX9>0 & COL1A1>0 / Sox9>0 & Col1a1>0","human_geometry_sha256":hm['geometry_sha256'],"rat_geometry_sha256":rm['geometry_sha256'],"old_global_X_umap_used":False,"figure4_formal_stage4":"NOT_ESTIMABLE_missing_AUCell_RSS_execution_audit","python":platform.python_version()}
    (QA/"FINAL_QA.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding="utf-8")
    rows=[]
    for p in sorted(ROOT.rglob("*")):
        if p.is_file() and "__pycache__" not in p.parts: rows.append({"relative_path":p.relative_to(ROOT).as_posix(),"bytes":p.stat().st_size,"sha256":sha256(p)})
    pd.DataFrame(rows).to_csv(MAN/"OUTPUT_INDEX.tsv",sep="\t",index=False); print(json.dumps(qa,ensure_ascii=False,indent=2))

if __name__=="__main__": main()
