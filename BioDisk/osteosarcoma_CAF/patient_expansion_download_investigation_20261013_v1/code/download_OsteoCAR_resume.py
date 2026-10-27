
from pathlib import Path
import os,json,urllib.request,urllib.error,concurrent.futures,threading,time,hashlib,datetime,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1');D=O/'downloads/OsteoCAR';files=json.loads((O/'resources/OSTEOCAR_DOWNLOAD_MANIFEST.json').read_text());chunk=32*1024*1024;lock=threading.Lock();handles={};tasks=[];results=[];failures=[];remaining={};file_records={}
def hashes(path):
 m=hashlib.md5();s=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):m.update(b);s.update(b)
 return m.hexdigest(),s.hexdigest()
prior=[]
for filename in ['OsteoCAR_chunks.jsonl','OsteoCAR_resume_chunks.jsonl']:
 f=O/'logs'/filename
 if f.exists():prior += [json.loads(l) for l in f.read_text().splitlines()]
completed={}
for f in files:
 final=D/f['name'];part=D/(f['name']+'.part')
 if final.exists():
  md5,sha=hashes(final);assert md5==f['computed_md5'];file_records[f['name']]=dict(path=str(final),bytes=final.stat().st_size,md5=md5,sha256=sha,official_md5_verified=True);continue
 fd=os.open(part,os.O_RDWR|os.O_CREAT,0o660);os.ftruncate(fd,f['size']);handles[f['name']]=fd
 for r in prior:
  if r['file']==f['name']:
   b=os.pread(fd,r['end']-r['start']+1,r['start']);assert hashlib.sha256(b).hexdigest()==r['sha256'];completed[f['name'],r['start'],r['end']]=True
 ft=[(f,s,min(s+chunk,f['size'])-1) for s in range(0,f['size'],chunk) if (f['name'],s,min(s+chunk,f['size'])-1) not in completed];tasks.extend(ft);remaining[f['name']]=len(ft)
print('RESUME_VALIDATED',len(completed),'chunks; remaining',len(tasks),flush=True)
def finalize(f):
 part=D/(f['name']+'.part');os.fsync(handles[f['name']]);md5,sha=hashes(part);assert md5==f['computed_md5'],f['name']+' MD5 mismatch';final=D/f['name'];part.rename(final);record=dict(path=str(final),bytes=final.stat().st_size,md5=md5,sha256=sha,official_md5_verified=True);file_records[f['name']]=record
 (O/'resources'/('OsteoCAR_'+f['name'].replace('.qs2','')+'_DOWNLOAD_COMPLETE.json')).write_text(json.dumps(record,indent=2));print('FILE_VERIFIED',f['name'],flush=True)
def fetch(task):
 f,start,end=task
 for attempt in range(1,7):
  try:
   op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:18797','http':'http://127.0.0.1:18797'}));req=urllib.request.Request(f['download_url'],headers={'Range':f'bytes={start}-{end}','Cache-Control':'no-cache'})
   with op.open(req,timeout=90) as resp:
    assert resp.status==206 and resp.headers.get('Content-Range')==f'bytes {start}-{end}/{f["size"]}'
    pos=start;h=hashlib.sha256()
    while True:
     b=resp.read(1048576)
     if not b:break
     assert pos+len(b)<=end+1;os.pwrite(handles[f['name']],b,pos);pos+=len(b);h.update(b)
    assert pos==end+1
   record=dict(file=f['name'],start=start,end=end,sha256=h.hexdigest(),attempt=attempt,completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
   with lock:
    with (O/'logs/OsteoCAR_resume_chunks.jsonl').open('a') as log:log.write(json.dumps(record)+'\n')
    remaining[f['name']]-=1;print('CHUNK_COMPLETE',f['name'],start,end,flush=True)
    if remaining[f['name']]==0:finalize(f)
   return None
  except Exception as e:
   error=type(e).__name__+': '+str(e)[:180]
   with lock:print('RETRY',f['name'],start,attempt,error,flush=True)
   if attempt<6:time.sleep(min(5*attempt,20))
 return dict(file=f['name'],start=start,end=end,error=error)
try:
 for f in files:
  if f['name'] in remaining and remaining[f['name']]==0:finalize(f)
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  futures=[pool.submit(fetch,t) for t in tasks]
  for future in concurrent.futures.as_completed(futures):
   result=future.result()
   if result:failures.append(result)
 for fd in handles.values():os.close(fd)
 handles={}
 if failures:
  (O/'resources/OSTEOCAR_REMAINING_FAILED_RANGES.json').write_text(json.dumps(failures,indent=2));raise RuntimeError(str(len(failures))+' ranges remain incomplete')
 assert len(file_records)==2
 (O/'resources/OSTEOCAR_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'status':'complete_verified','completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':list(file_records.values())},indent=2));(O/'logs/OsteoCAR_download.exit_code').write_text('0\n');print('ALL_OSTEOCAR_DOWNLOADS_VERIFIED',flush=True)
except BaseException as e:
 for fd in handles.values():os.close(fd)
 (O/'logs/OsteoCAR_download.exit_code').write_text('1\n');traceback.print_exc();raise
