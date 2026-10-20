from pathlib import Path
import hashlib,json,time,sys
from prepare_full_selected import ROOT,RNA,META,CLIN,SKILL
sys.stdout.reconfigure(encoding='utf8')
atlas=json.loads((ROOT/'audit/Figure1_GLOBAL_UMAP_QA.json').read_text(encoding='utf8'))
paths=[RNA,META,CLIN]+[Path(v['path']) for v in atlas['source_inputs'].values()]
records=[]
for p in paths:
    before=p.stat();hh=hashlib.sha256();read=0;last=time.time()
    print('Hash source',p.name,'bytes',before.st_size,flush=True)
    with p.open('rb') as f:
        for block in iter(lambda:f.read(16*1024*1024),b''):
            hh.update(block);read+=len(block)
            if time.time()-last>25:print(p.name,round(read/before.st_size*100,1),'percent hashed',flush=True);last=time.time()
    after=p.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    records.append({'path':str(p),'bytes':before.st_size,'mtime_ns':before.st_mtime_ns,'sha256':hh.hexdigest(),'unchanged_during_read':True})
(ROOT/'audit/PRIMARY_SOURCE_SHA256_MANIFEST.json').write_text(json.dumps({'status':'PASS_READ_ONLY_SOURCE_HASHING','sources':records},ensure_ascii=False,indent=2),encoding='utf8')
print('Full parent/source SHA-256 complete',flush=True)
