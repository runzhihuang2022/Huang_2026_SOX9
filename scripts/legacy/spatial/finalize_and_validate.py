from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import zipfile

import fitz
import numpy as np
import pandas as pd
from PIL import Image

sys.stdout.reconfigure(encoding='utf8')
ROOT=Path(__file__).resolve().parents[1]
PREVIOUS=Path('data/external/INPUT_9aef9219b6/20261009_Figure2AB_Figure1_matched_inputs_v02')
S1=Path('data/external/INPUT_f0bd868f67/Table_S1_Figure1_4_actual_scRNA_samples_v05.xlsx')
S1_EXPECTED='4303fe64c21b56bfe422fe6a259ba1d7b0f3b51b2224c9b37009007d7d36ffb6'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def validate():
    inp=json.loads((ROOT/'audit/INPUT_AND_MARKER_AUDIT.json').read_text(encoding='utf8'))
    qa=json.loads((ROOT/'audit/RENDER_AND_CANVAS_QA.json').read_text(encoding='utf8'))
    assert qa['palette']['double_positive']=='#E6B800' and qa['user_yellow_overlap_choice']
    assert all(q['outside_slot_pixel_changes']==0 for q in qa['composites'])
    assert qa['current_Figure1bc_copy_sha256']==qa['source_Figure1bc_sha256']
    protected=[]
    for r in inp['human_source_parts']+inp['rat_main_sections']:
        p=Path(r['source']);assert p.stat().st_size==r['bytes'] and p.stat().st_mtime_ns==r['mtime_ns']
        digest=sha(p);assert digest==r['sha256'];protected.append({'source':str(p),'sha256_unchanged':digest})
        name=r.get('data_key',r.get('sample'))
        with np.load(ROOT/'data'/f'{name}.npz') as d:
            fib=d['fib'];v=d['vals'][:,:2];ids=d['ids'];xy=d['xy']
            assert np.isfinite(v).all() and (v>=0).all() and np.isfinite(xy).all()
            assert len(ids)==len(pd.unique(ids))
            dual=fib&(v>0).all(axis=1);sox=fib&(v[:,0]>0);coll=fib&(v[:,1]>0)
            assert int(dual.sum())==r['double_positive_Fib']
            assert int(sox.sum())==r['SOX9_positive_Fib'] and int(coll.sum())==r['COL1A1_positive_Fib']
    previous=pd.read_csv(PREVIOUS/'tables/Figure3G_used_sections_metadata_counts.csv').set_index('sample')
    current=pd.read_csv(ROOT/'tables/Figure3G_main_sections_counts.csv').set_index('sample')
    assert set(previous.index)==set(current.index)
    for name,r in current.iterrows():
        assert r.double_positive_Fib==previous.loc[name,'Sox9_Col1a1_double_positive']
        assert r.target_in_registered_zoom==previous.loc[name,'target_in_registered_zoom']
    main=pd.read_csv(ROOT/'tables/Figure1E_main_sections_counts.csv').set_index('section')
    assert main.loc['2mph','double_positive_Fib']==149 and main.loc['6mph','double_positive_Fib']==152
    assert main.loc['2mph','target_in_registered_zoom']==24 and main.loc['6mph','target_in_registered_zoom']==24
    region=pd.read_csv(ROOT/'tables/Human_marker_by_physical_section_region.csv')
    assert set(region.section)=={'2 mph','6 mph','9 mph'}
    assert set(region.P_value)=={'NE'} and set(region.BH_FDR)=={'NE'}
    assert np.allclose(region.percent,region.dual_positive/region.n_source_Fib*100)
    # Verify typography from actual PDF text, not merely plotting rcParams.
    pdf_checks=[]
    for r in qa['figures']:
        base=ROOT/'figures'/r['name'];doc=fitz.open(base.with_suffix('.pdf'));assert len(doc)==1
        text=doc[0].get_text();spans=[]
        for block in doc[0].get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                spans.extend(line['spans'])
        sizes=[s['size'] for s in spans if s['text'].strip()]
        assert min(sizes)>=7.999 and max(sizes)<=12.001,(r['name'],sizes)
        assert all('Arial' in s['font'] for s in spans if s['text'].strip())
        fonts=doc[0].get_fonts(full=True)
        assert all(len(doc.extract_font(f[0])[3])>0 for f in fonts)
        actualmm=[doc[0].rect.width/72*25.4,doc[0].rect.height/72*25.4]
        assert np.allclose(actualmm,r['mm'],rtol=0,atol=.005)
        if 'unified_slot' in r['name']:
            assert 'PDGFRA' not in text and 'VIM' not in text and 'ACTA2' not in text
            assert ('SOX9' in text and 'COL1A1' in text) if 'Figure1E' in r['name'] else ('Sox9' in text and 'Col1a1' in text)
        with Image.open(base.with_suffix('.png')) as im:
            assert im.width>0 and im.height>0
            assert abs(im.info['dpi'][0]-600)<.1
        with Image.open(base.with_suffix('.tiff')) as im:
            assert tuple(float(d) for d in im.info['dpi'])==(600.,600.)
        pdf_checks.append({'name':r['name'],'actual_pdf_mm':actualmm,'actual_font_pt':[min(sizes),max(sizes)],'embedded_fonts':[f[3] for f in fonts],'gene_label_text':text})
        doc.close()
    assert sha(S1)==S1_EXPECTED,'Current scRNA Table S1 unexpectedly changed'
    for r in inp['reference_figures']:
        assert sha(r['source'])==r['sha256']
    for name in ['Figure1E_original_canvas_position_proof','Figure3_COLI_SOX9_updated_A4']:
        doc=fitz.open(ROOT/'figures'/f'{name}.pdf');assert len(doc)==1
        assert abs(doc[0].rect.width/72*25.4-210)<.005 and abs(doc[0].rect.height/72*25.4-297)<.005
        doc.close()
    final={'status':'PASS_NEW_DUAL_MARKER_PANELS_AND_UNCHANGED_SOURCES_WITH_EXPLICIT_LIMITS','user_confirmed_colors':{'SOX9':'red','COL I / COL1A1':'green','double_positive':'yellow'},'gate_recomputed_from_source_X':True,'human_main_dual_counts':{'2mph':149,'6mph_combined':152},'rat_counts_match_previous_dual_marker_analysis':True,'other_canvas_pixels_unchanged':True,'old_marker_names_removed_from_new_panel_text':True,'source_checks':protected,'PDF_typography_and_size':pdf_checks,'biological_replicate_test':'Not performed. P=NE; FDR=NE. Spatial observations not independent patients/animals.','Figure1_full_canvas_scope':'Position proof for e only; inherited b-c is old reference. Use approved current b-c separately.','Table_S1_scRNA_unchanged_sha256':S1_EXPECTED,'original_IF_pixels_not_generated_or_recolored':True,'source_resolution_caution':'Existing spatial observations preserved; cell segmentation has not been independently revalidated.'}
    (ROOT/'audit/FINAL_QA.json').write_text(json.dumps(final,ensure_ascii=False,indent=2),encoding='utf8')
    qa['status']='PASS_FINAL_DATA_PDF_AND_CANVAS_VALIDATION';(ROOT/'audit/RENDER_AND_CANVAS_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
    return final


class Links(HTMLParser):
    def __init__(self):super().__init__();self.urls=[]
    def handle_starttag(self,tag,attrs):
        self.urls.extend(v for k,v in attrs if k in ['src','href'])


def package():
    archive=ROOT/'outputs/SOX9_Figure1E_3G_COLI_SOX9_unified_v01.zip'
    checks=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and p.parent.name!='outputs' and '__pycache__' not in p.parts:
            checks.append({'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
    pd.DataFrame(checks).to_csv(ROOT/'audit/OUTPUT_CHECKSUMS.csv',index=False,encoding='utf-8-sig')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(ROOT.rglob('*')):
            if p.is_file() and p.parent.name!='outputs' and '__pycache__' not in p.parts:
                z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        archived=set(z.namelist())
        parser=Links();page=ROOT/'06_html/OPEN_ME_FIRST_Figure1E_3G_COLI_SOX9_unified.html';parser.feed(page.read_text(encoding='utf8'))
        for url in parser.urls:
            if '://' in url or url.startswith('#'):continue
            target=(page.parent/url).resolve();assert target.is_file(),url
            if target!=archive:assert target.relative_to(ROOT).as_posix() in archived,url
    print(json.dumps({'archive':str(archive),'archive_MB':round(archive.stat().st_size/1e6,2),'HTML_local_links_checked':len(parser.urls),'archive_CRC_verified':True},ensure_ascii=False),flush=True)


if __name__=='__main__':
    result=validate();print(result['status'],flush=True);package()
