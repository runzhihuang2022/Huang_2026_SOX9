#!/usr/bin/env python3
"""Opt-in source-H5AD spatial gate/region reconstruction, using approved readers.

No source RNA, barcode, individual-patient key, or full coordinates are exported.
Requires separately authorized external inputs; no network or SSH calls.
"""
from pathlib import Path
import argparse,ast,hashlib,json,sys,time
import numpy as np,pandas as pd,h5py
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',default='data/external/spatial')
    p.add_argument('--sources-json',help='Optional PRIVATE JSON mapping section IDs to authorized source paths; never commit it')
    p.add_argument('--out',required=True)
    p.add_argument('--maps',action='store_true',help='Also redraw full/ROI maps using the archived primitives')
    a=p.parse_args();out=Path(a.out).resolve()
    if out.exists():raise FileExistsError('Refusing to overwrite '+str(out))
    out.mkdir(parents=True);start=time.time()
    cfg=json.loads((ROOT/'config/spatial_sections.json').read_text(encoding='utf-8'))
    source_map=json.loads(Path(a.sources_json).read_text(encoding='utf-8-sig')) if a.sources_json else {}
    source=ROOT/'scripts/legacy/spatial/audit_and_prepare.py'
    tree=ast.parse(source.read_text(encoding='utf-8'))
    funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['sha','col','markers','rotate','read_source']]
    ns={'Path':Path,'np':np,'pd':pd,'h5py':h5py,'sparse':sparse,'hashlib':hashlib}
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(source),'exec'),ns)
    if a.maps:
        from matplotlib.font_manager import FontProperties,findfont
        findfont(FontProperties(family='Arial'),fallback_to_default=False)
        plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
        source=ROOT/'scripts/renderers/spatial_primitives.py'
        tree=ast.parse(source.read_text(encoding='utf-8'))
        exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),str(source),'exec'),ns)
        ns['Rectangle']=Rectangle
    rows=[];regions=[];checks=[]
    for s in cfg:
        ident=s['section_id'];path=Path(source_map.get(ident,Path(a.input_dir)/s['filename']))
        print('READ',ident,flush=True)
        ids,sub,xy,vals,rec=ns['read_source'](path,s['genes'],s['coordinate_key'])
        assert rec['sha256']==s['sha256'],'Source hash mismatch: '+ident
        fib=np.char.startswith(sub,'Fib_');target=fib&(vals[:,:2]>0).all(axis=1)
        for key in ['spatial_observations','source_Fib_observations','double_positive_Fib']:
            passed=rec[key]==s[key];checks.append({'section':ident,'check':key,'passed':passed})
            if not passed:raise AssertionError(ident+' '+key)
        row={k:rec[k] for k in ['spatial_observations','source_Fib_observations','SOX9_positive_Fib','COL1A1_positive_Fib',
                                'double_positive_Fib','SOX9_only_Fib','COL1A1_only_Fib','double_negative_Fib','double_positive_gt1','sha256']}
        row.update({'section_id':ident,'species':s['species'],'display_section':s['display_section'],
            'percent_dual_positive':100*int(target.sum())/int(fib.sum()),'P_value':'NE','BH_FDR':'NE'})
        if s['species']=='Human':
            xy=ns['rotate'](xy,s['rotation_degrees'])
            (x0,y0),(x1,y1)=np.asarray(s['anchor']['origin_line'],float)
            xb=np.full(len(xy),(x0+x1)/2) if np.isclose(y0,y1) else x0+(xy[:,1]-y0)*(x1-x0)/(y1-y0)
            region=np.where(xy[:,0]>xb,'Scarless side','Scar side')
            for r in ['Scarless side','Scar side']:
                m=fib&(region==r);regions.append({'section':s['display_section'],'source_part':ident,'region':r,
                    'n_source_Fib':int(m.sum()),'dual_positive':int((target&m).sum())})
        rows.append(row)
        if a.maps:
            roi=s['roi'];full0=ns['rectangular_limits'](xy)
            full=full0 if roi is None else (min(full0[0],roi[0]),max(full0[1],roi[1]),min(full0[2],roi[2]),max(full0[3],roi[3]))
            fig,ax=plt.subplots(1,2,figsize=(105/25.4,52.5/25.4));fig.subplots_adjust(left=.02,right=.98,bottom=.15,top=.85,wspace=.15)
            ns['plot_distribution'](ax[0],xy,target,full,'#E6B800',zoom=False)
            ns['add_scale_bar'](ax[0],full,.5)
            if roi is not None:
                ns['add_zoom_roi_box'](ax[0],roi)
                ns['plot_distribution'](ax[1],xy,target,roi,'#E6B800',zoom=True);ns['add_scale_bar'](ax[1],roi,.5)
            else:
                ax[1].axis('off');ax[1].text(.5,.5,'No source-part ROI\n(pooled-display ROI differs)',ha='center',va='center',fontsize=8)
            fig.suptitle(ident,fontsize=9);fig.text(.02,.02,'Dual positive: both RNA markers >0 in source Fib; other bins grey.',fontsize=8)
            for ext in ['pdf','svg','png']:fig.savefig(out/(ident+'.'+ext),dpi=600)
            plt.close(fig)
    pd.DataFrame(rows).to_csv(out/'source_level_spatial_counts.csv',index=False,encoding='utf-8-sig')
    bypart=pd.DataFrame(regions)
    expected=pd.read_csv(ROOT/'data/processed/spatial/Human_marker_by_source_part_region.csv').set_index(['section','source_part','region']).sort_index()
    actual=bypart.set_index(['section','source_part','region']).sort_index()
    for c in ['n_source_Fib','dual_positive']:
        ok=bool(np.array_equal(actual[c].to_numpy(),expected[c].to_numpy()));checks.append({'check':'human_regions_'+c,'passed':ok})
        if not ok:raise AssertionError('Human region counts mismatch')
    bypart.to_csv(out/'human_source_regions.csv',index=False,encoding='utf-8-sig')
    report={'status':'executed_verified','scope':'Source-H5AD marker gate and registered human regions; map panels optional; no biological-replicate statistical claim',
            'elapsed_seconds':time.time()-start,'checks':checks,'section_count':len(rows),'maps_rerendered':a.maps}
    (out/'source_spatial_run_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'status':report['status'],'sections':len(rows),'checks':len(checks),'elapsed_seconds':report['elapsed_seconds']}))

if __name__=='__main__':main()
