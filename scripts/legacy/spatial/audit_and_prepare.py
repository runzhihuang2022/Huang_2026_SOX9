from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from PIL import Image

sys.stdout.reconfigure(encoding='utf8')
ROOT = Path(__file__).resolve().parents[1]
OLD = Path('data/external/INPUT_03d0736bfe/13_sequence_analysis_SOX9_COL1A1_cutoff0_20260917_v05')
HROIS = Path('data/external/INPUT_577c4659fd/FIGURE_INDEX.tsv')
HPARENT = Path('data/external/INPUT_d05dd65597/02_canonical')
RROIS = Path('data/external/INPUT_d28991e41f/FIGURE_INDEX.tsv')
RAT = Path('data/external/INPUT_01762fefba/spatial_h5ad_compatible')
SOURCE = Path('data/external/INPUT_6ab40abaef/SOX9投稿20261008')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(b)
    return digest.hexdigest()


def col(g):
    if isinstance(g, h5py.Dataset):
        return g.asstr()[:] if g.dtype.kind in 'OSU' else g[:]
    if 'values' in g:
        a = col(g['values']).astype(object)
        if 'mask' in g:
            a[g['mask'][:]] = None
        return a
    c = g['codes'][:]
    a = col(g['categories'])[np.maximum(c, 0)].astype(object)
    a[c < 0] = None
    return a


def markers(h, genes):
    vg = h['var']; names = col(vg[vg.attrs['_index']])
    positions = [np.flatnonzero(names == gene) for gene in genes]
    assert all(len(p) == 1 for p in positions), genes
    positions = np.array([int(p[0]) for p in positions])
    x = h['X']; shape = tuple(x.attrs['shape'])
    out = np.zeros((shape[0], len(genes)), dtype=np.float64)
    if x.attrs['encoding-type'] == 'csr_matrix':
        ptr = x['indptr'][:]
        for a in range(0, shape[0], 2048):
            b = min(a + 2048, shape[0]); p0 = int(ptr[a]); p1 = int(ptr[b])
            block = sparse.csr_matrix((x['data'][p0:p1], x['indices'][p0:p1], ptr[a:b+1] - p0), shape=(b-a, shape[1]))
            out[a:b] = block[:, positions].toarray()
    else:
        assert x.attrs['encoding-type'] == 'csc_matrix'
        ptr = x['indptr'][:]
        for k, j in enumerate(positions):
            a, b = int(ptr[j]), int(ptr[j+1])
            np.add.at(out[:, k], x['indices'][a:b], x['data'][a:b])
    assert np.isfinite(out).all() and (out >= 0).all()
    return out


def rotate(xy, angle, center=None):
    t = np.deg2rad(angle); r = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    c = np.zeros(2) if center is None else np.asarray(center)
    return (xy - c) @ r.T + c


def read_source(path, genes, key):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    digest = sha(path)
    with h5py.File(path, 'r') as h:
        og = h['obs']; ids = col(og[og.attrs['_index']]).astype(str)
        assert pd.Index(ids).is_unique
        sub = col(og['sub_labels']).astype(str)
        xy = h['obsm'][key][:].astype(float)
        vals = markers(h, genes)
        unit_evidence = {k: str(h['uns'][k][()]) for k in ['bin_size', 'cellbin', 'resolution'] if k in h['uns'] and isinstance(h['uns'][k], h5py.Dataset)}
        rec = {'source': str(path), 'sha256': digest, 'bytes': before[0], 'mtime_ns': before[1], 'matrix_shape': [int(n) for n in h['X'].attrs['shape']], 'matrix': 'X, nonnegative normalized RNA as stored; no re-normalization', 'obsm_keys': list(h['obsm']), 'coordinate_key': key, 'spatial_unit_metadata': unit_evidence, 'log1p_metadata_present': 'log1p' in h['uns']}
    assert len(ids) == len(xy) == len(vals) and np.isfinite(xy).all()
    assert before == (path.stat().st_size, path.stat().st_mtime_ns)
    fib = np.char.startswith(sub, 'Fib_')
    sox, coll = vals[:, 0] > 0, vals[:, 1] > 0
    rec.update({'spatial_observations': len(ids), 'source_Fib_observations': int(fib.sum()), 'SOX9_positive_Fib': int((fib & sox).sum()), 'COL1A1_positive_Fib': int((fib & coll).sum()), 'double_positive_Fib': int((fib & sox & coll).sum()), 'SOX9_only_Fib': int((fib & sox & ~coll).sum()), 'COL1A1_only_Fib': int((fib & coll & ~sox).sum()), 'double_negative_Fib': int((fib & ~sox & ~coll).sum()), 'double_positive_gt1': int((fib & (vals[:, :2] > 1).all(axis=1)).sum()), 'minimum_positive_SOX9': float(vals[:, 0][sox].min()) if sox.any() else None, 'minimum_positive_COL1A1': float(vals[:, 1][coll].min()) if coll.any() else None})
    assert rec['double_positive_Fib'] + rec['SOX9_only_Fib'] + rec['COL1A1_only_Fib'] + rec['double_negative_Fib'] == rec['source_Fib_observations']
    return ids, sub, xy, vals, rec


def save_data(name, ids, sub, xy, vals, fib, **kwargs):
    np.savez_compressed(ROOT/'data'/f'{name}.npz', ids=ids, sub=sub, xy=xy, vals=vals, fib=fib, **kwargs)


def main():
    for d in ['audit', 'data', 'tables', 'figures', '06_html', 'outputs']:
        (ROOT/d).mkdir(parents=True, exist_ok=True)
    human_manifest = json.loads((OLD/'00_manifest/Figure1E_human_spatial_manifest_v05.json').read_text(encoding='utf8'))
    rat_manifest = json.loads((OLD/'00_manifest/Figure3G_rat_spatial_manifest_v05.json').read_text(encoding='utf8'))
    hrois = pd.read_csv(HROIS, sep='\t').set_index('section')
    rrois = pd.read_csv(RROIS, sep='\t').set_index('sample')
    human=[]; arrays={}
    for entry in human_manifest['inputs']:
        path = Path(entry['input']); name = Path(entry['anchor']).name.replace('_spatial_anchors.json', '')
        ids, sub, raw, vals, rec = read_source(path, ['SOX9', 'COL1A1', 'PDGFRA', 'VIM'], 'spatial')
        assert rec['sha256'] == entry['sha256'], name
        anchor = json.loads(Path(entry['anchor']).read_text(encoding='utf8'))
        assert sha(entry['anchor']) == entry['anchor_sha256']
        fib = np.char.startswith(sub, 'Fib_'); angle = float(anchor.get('rotation_applied', 0))
        display = rotate(raw, angle); line = np.asarray(anchor['origin_line'], float)
        (x0, y0), (x1, y1) = line
        xbound = np.full(len(raw), (x0+x1)/2) if np.isclose(y0, y1) else x0 + (display[:,1]-y0)*(x1-x0)/(y1-y0)
        assert anchor['analysis_mode'] == 'scar_to_scarless' and anchor['target_direction'] == 'right_is_scarless'
        region = np.where(display[:,0] > xbound, 'Scarless side', 'Scar side')
        legacy = fib & (vals[:,0] > 1)
        rec.update({'species':'Human', 'section':entry['section'], 'data_key': name, 'genes':['SOX9', 'COL1A1'], 'rotation_degrees':angle, 'anchor':entry['anchor'], 'anchor_sha256':entry['anchor_sha256'], 'anchor_rotation_reference':'rotation about origin; coordinates as recorded in origin_line', 'legacy_SOX9_gt1_gate':int(legacy.sum()), 'legacy_3markers_all_gt0':int((fib & (vals[:,[0,2,3]]>0).all(axis=1)).sum()), 'current_main_usage':name in ['2mph','6mph_p1','6mph_p2'], 'independent_donor_id':'Not recorded', 'statistical_status':'Spatial observations not independent replicates; 6 mph parts from same named source section'})
        arrays[name] = (ids, sub, raw, vals, fib, display, region, line, rec)
        save_data(name, ids, sub, display, vals, fib, region=region, boundary=line)
        human.append(rec)
        print('AUDIT HUMAN', name, rec['spatial_observations'], rec['double_positive_Fib'], flush=True)
    # Reuse EXACT registered display geometry from the old human main-figure atlas.
    # This atlas concatenates the 6 mph parts without new translation or re-selection.
    main_human=[]
    for section, names in [('2mph',['2mph']), ('6mph',['6mph_p1','6mph_p2'])]:
        parent = HPARENT/section/f'{section}.sox9_fib_dual_metadata.h5ad'
        with h5py.File(parent,'r') as h:
            parent_ids=col(h['obs'][h['obs'].attrs['_index']]).astype(str)
            parent_raw=h['obsm/spatial_source'][:].astype(float)
        wanted=np.concatenate([np.char.add(name+'::', arrays[name][0]) for name in names])
        loc=pd.Index(wanted).get_indexer(parent_ids); assert (loc>=0).all() and len(np.unique(loc))==len(wanted)
        raw=np.vstack([arrays[name][2] for name in names])[loc]
        assert np.allclose(raw,parent_raw,rtol=0,atol=0),section+' source coordinate mismatch'
        vals=np.vstack([arrays[name][3] for name in names])[loc]
        sub=np.concatenate([arrays[name][1] for name in names])[loc]
        fib=np.char.startswith(sub,'Fib_')
        angle=45 if section=='2mph' else 0
        # The v04 plotting adapter used the reference origin-rotation on spatial_source.
        display=rotate(parent_raw,angle)
        r=hrois.loc[section]; roi=np.array([r[k] for k in ['roi_xmin','roi_xmax','roi_ymin','roi_ymax']],float)
        inside=(display[:,0]>=roi[0])&(display[:,0]<=roi[1])&(display[:,1]>=roi[2])&(display[:,1]<=roi[3])
        target=fib&(vals[:,:2]>0).all(axis=1)
        rec={'species':'Human','section':section, 'source_parts':names,'spatial_observations':len(parent_ids),'source_Fib_observations':int(fib.sum()), 'double_positive_Fib':int(target.sum()), 'target_in_registered_zoom':int((target&inside).sum()), 'roi':roi.tolist(), 'coordinate_key':'spatial_source, origin-rotation', 'rotation_degrees':angle, 'source_geometry_parent':str(parent), 'source_geometry_sha256':sha(parent), 'um_per_coordinate_unit':0.5,'zoom_selection':'Reuse old SOX9-only registered ROI, not re-tuned to dual marker'}
        save_data('main_'+section,parent_ids,sub,display,vals,fib,roi=roi)
        main_human.append(rec); print('MAIN HUMAN',section,rec,flush=True)
    rats=[]
    sections=[('Adult','1 dpw','Adult_1d_C01830C1D1'),('Adult','5 dpw','Adult_5d1_C01830E6F6'),('Adult','19 dpw','Adult_19d_B03424F4'),('E16.5','6 h','E16.5_6H_C02938E1'),('E16.5','3 dpw','E16.5_3D_D03056D5'),('E16.5','6 dpw','E16.5_6D_C02938F4')]
    expected={x['sample']:x['sha256'] for x in rat_manifest['inputs']}
    for stage,time,sample in sections:
        path=RAT/f'{sample}_Annotated_Spatial.h5ad'
        ids,sub,xy,vals,rec=read_source(path,['Sox9','Col1a1'],'X_rot_spatial')
        assert rec['sha256']==expected[sample]
        r=rrois.loc[sample];assert r.coordinate_key=='X_rot_spatial'
        roi=np.array([r[k] for k in ['roi_xmin','roi_xmax','roi_ymin','roi_ymax']],float)
        fib=np.char.startswith(sub,'Fib_');target=fib&(vals>0).all(axis=1)
        inside=(xy[:,0]>=roi[0])&(xy[:,0]<=roi[1])&(xy[:,1]>=roi[2])&(xy[:,1]<=roi[3])
        rec.update({'species':'Rat','stage':stage,'wound_time':time,'sample':sample, 'genes':['Sox9','Col1a1'],'target_in_registered_zoom':int((target&inside).sum()),'roi':roi.tolist(),'additional_rotation':0,'um_per_coordinate_unit':0.5,'independent_animal_id':'Not recorded','statistical_status':'One registered section per displayed time; descriptive only'})
        rats.append(rec);save_data(sample,ids,sub,xy,vals,fib,roi=roi)
        print('AUDIT RAT',sample,rec['double_positive_Fib'],flush=True)
    pd.DataFrame(human).to_csv(ROOT/'tables/Human_four_source_parts_marker_audit.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(main_human).to_csv(ROOT/'tables/Figure1E_main_sections_counts.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(rats).to_csv(ROOT/'tables/Figure3G_main_sections_counts.csv',index=False,encoding='utf-8-sig')
    sources=[]
    for number in [1,3]:
        path=SOURCE/f'Figure {number}.tif'
        with Image.open(path) as im:
            rec={'figure':number,'source':str(path),'pixel_size':list(im.size),'dpi':[float(x) for x in im.info.get('dpi',[])],'sha256':sha(path)}
        sources.append(rec)
    def safe(x):
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x))
    payload={'gate_human':'source Fib_ AND SOX9>0 AND COL1A1>0','gate_rat':'source Fib_ AND Sox9>0 AND Col1a1>0','human_source_parts':human,'human_main_sections':main_human,'rat_main_sections':rats,'reference_figures':sources,'source_markers_not_protein_measurements':True,'IF_colors_in_reference':{'SOX9':'red','COL I':'green','DAPI':'blue'},'ST_requested_dual_population_color':'yellow, as explicitly selected by user; not synthetic IF merge','no_original_data_modified':True}
    (ROOT/'audit/INPUT_AND_MARKER_AUDIT.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,default=safe),encoding='utf8')
    print('INPUT AUDIT COMPLETE',flush=True)


if __name__=='__main__':main()
