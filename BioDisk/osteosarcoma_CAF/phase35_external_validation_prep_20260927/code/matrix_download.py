
from pathlib import Path
import requests,json,hashlib,datetime,threading,traceback,time,gzip
import pandas as pd
import h5py,numpy as np
from concurrent.futures import ThreadPoolExecutor,as_completed
O=Path(__file__).resolve().parents[1];jobs=pd.read_csv(O/'tables/matrix_download_manifest.tsv',sep='\t').fillna('').to_dict('records');lock=threading.Lock();done=[];active={}
def now():return datetime.datetime.now().astimezone().isoformat()
def state():
 with lock:
  f=O/'MATRIX_DOWNLOAD_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps({'at':now(),'scheduled':len(jobs),'finished':len(done),'active':dict(active),'completed':sum(x['status']=='complete' for x in done),'failed':sum(x['status']!='complete' for x in done)},indent=2));t.replace(f)
def fetch(j):
 f=Path(j['path']);f.parent.mkdir(parents=True,exist_ok=True);part=f.with_name(f.name+'.part');s=requests.Session();s.trust_env=False
 with lock:active[f.name]='downloading'
 state();error=''
 for attempt in range(3):
  try:
   if not f.exists():
    with s.get(j['url'],timeout=(15,120),stream=True) as r:
     r.raise_for_status();expected=r.headers.get('Content-Length');size=0
     with part.open('wb') as h:
      for b in r.iter_content(1048576):h.write(b);size+=len(b)
     if expected:assert size==int(expected)
     assert size>0;part.replace(f)
   h=hashlib.sha256();
   with f.open('rb') as r:
    for b in iter(lambda:r.read(8388608),b''):h.update(b)
   qa={}
   if j['role']=='matrix':
    with h5py.File(f,'r') as r:
     m=r['matrix'];sh=m['shape'][:];bc=m['barcodes'].asstr()[:];genes=m['features/name'].asstr()[:];v=m['data'][:]
     assert len(bc)==sh[1] and len(set(bc))==len(bc) and len(genes)==sh[0]
     assert np.all(np.isfinite(v)) and np.all(v>=0) and np.all(v==np.floor(v))
     assert m['indptr'].shape[0]==sh[1]+1 and m['indptr'][-1]==len(v)
     qa=dict(features=int(sh[0]),spots=int(sh[1]),nnz=len(v),counts_valid=True,measured_candidate_genes=[x for x in ['PDGFRB','MCAM','GPC3','THY1','FAP','JAG1','CSPG4','GPRC5A'] if x in set(genes)])
   elif f.suffix=='.gz':
    with gzip.open(f,'rb') as r:
     while r.read(1048576):pass
   with lock:active.pop(f.name,None)
   return dict(**j,status='complete',bytes=f.stat().st_size,sha256=h.hexdigest(),qa=qa)
  except Exception as e:error=repr(e)
 with lock:active.pop(f.name,None)
 return dict(**j,status='failed',error=error)
state()
with ThreadPoolExecutor(max_workers=4) as pool:
 for f in as_completed([pool.submit(fetch,j) for j in jobs]):
  with lock:done.append(f.result())
  state()
(O/'delivery/MATRIX_DOWNLOAD_OUTCOMES.json').write_text(json.dumps(done,indent=2))
pd.DataFrame([{**{k:v for k,v in x.items() if k!='qa'},'qa':json.dumps(x.get('qa',{}))} for x in done]).to_csv(O/'tables/matrix_download_integrity_manifest.tsv',sep='\t',index=False)
q=dict(at=now(),complete=sum(x['status']=='complete' for x in done),failed=sum(x['status']!='complete' for x in done),matrix_files=sum(x['role']=='matrix' and x['status']=='complete' for x in done),validation_analysis_performed=False)
(O/'delivery/MATRIX_DOWNLOAD_QA.json').write_text(json.dumps(q,indent=2));(O/'matrix_download.exit_code').write_text('0' if q['failed']==0 else '1');(O/'MATRIX_DOWNLOAD_FINISHED').write_text(now());print(json.dumps(q))
