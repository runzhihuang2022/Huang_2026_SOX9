#!/usr/bin/env python3
"""Read-only public-package checks. Not proof of consent, licenses, or science."""
from pathlib import Path
import ast, csv, hashlib, json, re, sys

ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def main():
    errors=[];files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                     and not any(x.startswith(('my_','ci_','full_')) for x in p.relative_to(ROOT).parts[:-1])]
    text_suffix={'.py','.md','.json','.csv','.tsv','.txt','.yml','.yaml','.svg','.html'}
    forbidden=[('Windows drive path',r'\b[A-Za-z]:[\\/]'),('private server path',r'/media/iq/'),
               ('private network',r'\b192\.168\.\d+\.\d+\b'),('source patient ID',r'\bpatient\d+\b'),
               ('private key',r'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----'),
               ('GitHub token',r'\b(?:ghp_|github_pat_)[A-Za-z0-9_]{20,}')]
    for p in files:
        rel=p.relative_to(ROOT).as_posix()
        if p.stat().st_size>=50*1024*1024:errors.append({'file':rel,'issue':'File >= 50 MiB; do not commit'})
        if p.suffix.lower() in {'.h5ad','.rds','.fastq','.bam','.loom','.pem','.key','.zip'}:errors.append({'file':rel,'issue':'Raw/credential/archive file included'})
        if p.suffix=='.py':
            try:ast.parse(p.read_text(encoding='utf-8-sig'))
            except SyntaxError as e:errors.append({'file':rel,'issue':'Python syntax error: '+str(e)})
        if p.suffix.lower() in text_suffix:
            text=p.read_text(encoding='utf-8-sig',errors='replace')
            # Validator contains intentional private-path patterns, not actual data.
            if rel=='scripts/validate_release.py':continue
            for label,pattern in forbidden:
                m=re.search(pattern,text)
                if m:errors.append({'file':rel,'issue':label,'line':text[:m.start()].count('\n')+1})
    manifest=ROOT/'provenance/release_checksums.csv'
    if manifest.exists() and '--skip-checksums' not in sys.argv:
        for r in rows(manifest):
            p=ROOT/r['path']
            if not p.exists() or sha(p)!=r['sha256']:errors.append({'file':r['path'],'issue':'Checksum mismatch'})
    s1=rows(ROOT/'data/metadata/Table_S1_current_public_summary.csv')
    if len(s1)!=104:errors.append({'issue':'Table S1 not 104 actual libraries'})
    if any('PBMC' in r['True tissue'].upper() for r in s1):errors.append({'issue':'Unexpected PBMC specimen'})
    whitelist=rows(ROOT/'config/skin_19_record_whitelist.csv')
    if len(whitelist)!=19:errors.append({'issue':'Skin whitelist not 19 records'})
    fib=sum(int(float(r['Figure2A/B Fib cells'])) for r in s1)
    if fib!=20270:errors.append({'issue':'Figure2 Fib total not 20270','actual':fib})
    code=rows(ROOT/'provenance/code_manifest.csv')
    if len(code)<20:errors.append({'issue':'Source-code inventory incomplete'})
    result={'status':'PASS' if not errors else 'FAIL','files_checked':len(files),'TableS1_rows':len(s1),
            'skin_records':len(whitelist),'Figure2_Fib_cells':fib,'errors':errors,
            'scope':'Paths, identifiers, sizes, syntax, pinned file hashes and sample counts; not consent/license clearance'}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0 if not errors else 1

if __name__=='__main__':sys.exit(main())
