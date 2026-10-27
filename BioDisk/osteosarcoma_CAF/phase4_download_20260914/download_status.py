"""SYNTAX_CHECK: python -m py_compile download_status.py
RUN_COMMAND: python download_status.py --root /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
EXPECTED_SUCCESS_MARKER: DOWNLOAD_STATUS_RECORDED
"""
import argparse,json,os,requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
manifest=json.loads((root/'download_manifest.json').read_text())
rows=[]
for s in manifest:
 for url in s['urls']:
  f=root/'downloads'/s['dataset']/url.rsplit('/',1)[-1];meta=f.with_name(f.name+'.verified.json')
  rows.append({'dataset':s['dataset'],'sample':s['sample'],'name':f.name,'url':url,'verified':f.exists() and meta.exists(),'complete_bytes':f.stat().st_size if f.exists() else 0,'partial_bytes':f.with_name(f.name+'.geo.part').stat().st_size if f.with_name(f.name+'.geo.part').exists() else 0})
def size(r):
 if r['verified']:r['expected_bytes']=r['complete_bytes'];return
 try:
  z=requests.head(r['url'],timeout=15,allow_redirects=True);z.raise_for_status();r['expected_bytes']=int(z.headers['Content-Length'])
 except Exception as e:r['expected_bytes']=None;r['size_error']=type(e).__name__
with ThreadPoolExecutor(max_workers=5) as pool:list(pool.map(size,rows))
pid=int((root/'python.pid').read_text())
result={'files_total':len(rows),'files_verified':sum(r['verified'] for r in rows),'verified_bytes':sum(r['complete_bytes'] for r in rows if r['verified']),'expected_bytes_known':sum(r['expected_bytes'] or 0 for r in rows),'unknown_sizes':sum(r['expected_bytes'] is None for r in rows),'python_pid':pid,'process_exists':Path('/proc',str(pid)).exists(),'remaining':[r for r in rows if not r['verified']]}
(root/'download_status.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print('DOWNLOAD_STATUS_RECORDED')
