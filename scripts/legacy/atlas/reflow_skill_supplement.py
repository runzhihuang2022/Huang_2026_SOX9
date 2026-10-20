from pathlib import Path
import ast, importlib.util, json, sys
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from prepare_full_selected import ROOT,SKILL,sha
sys.stdout.reconfigure(encoding='utf8')
p=SKILL/'scripts/run_figure2_ae_universal.py';source=p.read_text(encoding='utf8');tree=ast.parse(source)
main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='palette' for t in n.targets))
end=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='base' for t in n.targets))
run=ROOT/'analysis/human_universal';output=ROOT/'analysis/human_universal_presentation'
for folder in ['04_figures','05_tables']:(output/folder).mkdir(parents=True,exist_ok=True)
cfg=json.loads((ROOT/'00_manifest/figure2_config.json').read_text(encoding='utf8'))
manifest=json.loads((run/'00_manifest/Figure2_A-E_manifest.json').read_text(encoding='utf8'))
g=pd.read_csv(manifest['geometry_path'],sep='\t',float_precision='round_trip');clusters=manifest['clusters']
marker=pd.read_csv(run/'05_tables/Figure2_marker_allowlist_summary.tsv',sep='\t');ledger=pd.read_csv(run/'05_tables/Figure2_selected_marker_ledger.tsv',sep='\t')
mandatory=cfg['markers']['mandatory'];displayed=list(dict.fromkeys(mandatory+ledger.gene.tolist()))
spec=importlib.util.spec_from_file_location('pristine',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
module.setup_style(cfg['style'])
env=dict(module.__dict__);env.update({'cfg':cfg,'output':output,'cluster_order':clusters,'labels':g.cluster.to_numpy(),'coordinates':g[['UMAP1','UMAP2']].to_numpy(),'marker':marker,'displayed_genes':displayed,'mandatory':mandatory,'sample_summary':pd.read_csv(run/'05_tables/Figure2_sample_cluster_summary.tsv',sep='\t'),'seed':manifest['seed']})
exec(compile(ast.Module(body=main.body[start:end],type_ignores=[]),str(p),'exec'),env)
fig=env['fig'];axis=env['dot_axis']
for tick in axis.get_xticklabels():tick.set_rotation(90);tick.set_horizontalalignment('center')
# Only the overlapping x labels are reflowed; exact analytical/rendering bodies are executed unchanged.
base=output/'04_figures/Figure2_selected19_full_A-E_skill_layout_reflow'
for suffix in ['pdf','png','svg']:fig.savefig(base.with_suffix('.'+suffix),dpi=600)
plt.close(fig)
(ROOT/'audit/Supplement_AE_layout_adapter.json').write_text(json.dumps({'pristine_skill_sha256':sha(p),'plotting_AST_body_unchanged':True,'analytical_tables_unchanged':True,'only_layout_change':'Rotate 17-cluster dotplot x labels 90 degrees for readability','unrequested_E_status':'NOT_ESTIMABLE: no second preregistered contrast','frozen_geometry_sha256':manifest['geometry_sha256']},indent=2),encoding='utf8')
print('Skill supplementary A-E presentation reflowed; analytical source preserved',flush=True)
