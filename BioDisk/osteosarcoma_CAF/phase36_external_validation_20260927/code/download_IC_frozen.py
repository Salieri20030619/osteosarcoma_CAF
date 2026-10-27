
from pathlib import Path
import requests,json,gzip,hashlib,datetime
O=Path(__file__).resolve().parents[1]
out=[]
for name in ['GSE293065_ica.cell.embeddings.all.csv.gz','GSE293065_ica.feature.loadings.all.csv.gz']:
 f=O/'downloads'/name;part=f.with_suffix(f.suffix+'.part');url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293065/suppl/'+name
 s=requests.Session();s.trust_env=False
 for attempt in range(3):
  try:
   if not f.exists():
    with s.get(url,stream=True,timeout=(15,120)) as r:
     r.raise_for_status();expected=int(r.headers.get('Content-Length',0));size=0
     with part.open('wb') as h:
      for b in r.iter_content(1048576):
       h.write(b);size+=len(b)
       (O/'IC_DOWNLOAD_STATUS.json').write_text(json.dumps({'at':datetime.datetime.now().astimezone().isoformat(),'file':name,'bytes':size,'expected':expected,'completed':len(out)}))
     if expected:assert size==expected
     part.replace(f)
   h=hashlib.sha256(f.read_bytes())
   with gzip.open(f,'rt') as r:
    header=r.readline();n=sum(1 for _ in r)
   out.append(dict(file=str(f),source=url,bytes=f.stat().st_size,sha256=h.hexdigest(),rows=n,header=header[:1000],status='complete'));break
  except Exception as e:
   if attempt==2:out.append(dict(file=str(f),source=url,status='failed',error=repr(e)))
(O/'delivery/IC_DOWNLOAD_QA.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));raise SystemExit(0 if all(x['status']=='complete' for x in out) else 1)
