from pathlib import Path
import json,sys,zipfile,tempfile,hashlib
from html.parser import HTMLParser
from prepare_full_selected import ROOT,sha
from build_delivery_and_QA import local_files,Links
sys.stdout.reconfigure(encoding='utf8')
entry=ROOT/'06_html/OPEN_ME_FIRST_SOX9_selected19_full_Figure1_2_complete_v01.html'
archive=ROOT/'outputs/SOX9_Figure1BC_Figure2AB_selected19_fullcells_v01.zip'
assert (ROOT/'audit/PRIMARY_SOURCE_SHA256_MANIFEST.json').exists()
manifest=[]
for p in local_files(ROOT):
    if p==archive or p.name=='DELIVERY_FILE_MANIFEST.json' or p.suffix in ['.h5ad','.npy','.pyc'] or '__pycache__' in p.parts:continue
    manifest.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
(ROOT/'audit/DELIVERY_FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
files=[ROOT/x['relative_path'] for x in manifest]+[ROOT/'audit/DELIVERY_FILE_MANIFEST.json']
self_link='<a href="../outputs/SOX9_Figure1BC_Figure2AB_selected19_fullcells_v01.zip">离线图册与代码压缩包</a>'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:
        relative=p.relative_to(ROOT).as_posix()
        if p==entry:
            # A ZIP cannot include itself; extracted gallery has no broken self-download link.
            text=p.read_text(encoding='utf8').replace(self_link,'离线图册已解压')
            z.writestr(relative,text)
        else:z.write(p,relative)
with tempfile.TemporaryDirectory(prefix='SOX9_selected19_archive_QA_') as tmp:
    dest=Path(tmp)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        z.extractall(dest)
    page=dest/entry.relative_to(ROOT);parser=Links();parser.feed(page.read_text(encoding='utf8'))
    missing=[link for link in parser.links if not (page.parent/link).resolve().exists()]
    assert not missing,missing
    # All gallery links resolve within the fresh extraction without F-drive access.
    assert all((page.parent/link).resolve().is_relative_to(dest.resolve()) for link in parser.links)
    archive_links=len(parser.links)
live=Links();live.feed(entry.read_text(encoding='utf8'))
assert all((entry.parent/link).resolve().exists() for link in live.links)
qa={'status':'PASS_OFFLINE_EXTRACTED_GALLERY_LINKS_AND_ZIP_CRC','archive':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),'files':len(files),'extracted_gallery_links':archive_links,'missing_extracted_links':0,'live_gallery_links':len(live.links),'extracted_gallery_self_download_link_removed':True,'large_H5AD_and_PCA_objects_omitted_with_paths_hashes_and_staging_instructions':True}
(ROOT/'audit/OFFLINE_ARCHIVE_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(qa,ensure_ascii=False,indent=2),flush=True)
