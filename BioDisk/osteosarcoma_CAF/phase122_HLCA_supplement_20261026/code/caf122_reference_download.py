from pathlib import Path
import requests,concurrent.futures,time,json,os,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');url='https://datasets.cellxgene.cziscience.com/5f863718-02a6-46cb-83e7-52be988b915b.h5ad';total=21860171922;step=16*1024*1024
dest=O/'reference/HLCA_full_expression.h5ad.partial';journal=O/'logs/reference_download_chunks.jsonl';done={}
if journal.exists():
 for line in journal.read_text().splitlines():
  r=json.loads(line);done[r['start']]=r
if not dest.exists():
 assert not done
 with dest.open('wb') as f:f.truncate(total)
assert dest.stat().st_size==total
def fetch(start):
 end=min(start+step,total)-1
 for attempt in range(3):
  try:
   r=requests.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},timeout=(15,90));r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{total}';assert len(r.content)==end-start+1
   with dest.open('r+b',buffering=0) as f:f.seek(start);f.write(r.content)
   return dict(start=start,end=end,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
  except Exception:
   if attempt==2:raise
   time.sleep(2*(attempt+1))
pending=[s for s in range(0,total,step) if s not in done];t=time.time();written=0
with journal.open('a',encoding='utf-8') as f,concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for rec in ex.map(fetch,pending):
  f.write(json.dumps(rec)+'\n');f.flush();done[rec['start']]=rec;written+=rec['bytes']
  if len(done)%4==0:print('DOWNLOAD',len(done),'of',(total+step-1)//step,'MBps',round(written/max(time.time()-t,1)/1e6,2),flush=True)
assert len(done)==len(range(0,total,step))
hasher=hashlib.sha256()
with dest.open('rb') as f:
 for start in range(0,total,step):
  data=f.read(min(step,total-start));assert hashlib.sha256(data).hexdigest()==done[start]['sha256'];hasher.update(data)
final=dest.with_suffix('');dest.rename(final)
(O/'resources/REFERENCE_FULL_DOWNLOAD_COMPLETE.json').write_text(json.dumps(dict(status='PASS',url=url,bytes=total,sha256=hasher.hexdigest(),chunks=len(done),all_saved_chunks_rehashed=True),indent=2),encoding='utf-8');print('COMPLETE',flush=True)
