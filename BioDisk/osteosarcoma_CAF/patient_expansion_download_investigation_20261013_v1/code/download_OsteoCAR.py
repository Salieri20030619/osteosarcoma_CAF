
from pathlib import Path
import os,json,urllib.request,concurrent.futures,threading,time,hashlib,datetime,shutil,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1');D=O/'downloads/OsteoCAR';files=json.loads((O/'resources/OSTEOCAR_DOWNLOAD_MANIFEST.json').read_text());chunk=32*1024*1024
lock=threading.Lock();progress=O/'logs/OsteoCAR_chunks.jsonl';opener=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:18797','https':'http://127.0.0.1:18797'}));handles={};tasks=[]
def digest_file(f):
 m=hashlib.md5();s=hashlib.sha256()
 with f.open('rb') as h:
  for b in iter(lambda:h.read(8*1024*1024),b''):m.update(b);s.update(b)
 return m.hexdigest(),s.hexdigest()
try:
 assert shutil.disk_usage(O).free>sum(f['size'] for f in files)+5_000_000_000,'Insufficient disk capacity'
 for f in files:
  path=D/(f['name']+'.part');fd=os.open(path,os.O_RDWR|os.O_CREAT,0o660);os.ftruncate(fd,f['size']);handles[f['name']]=fd
  for start in range(0,f['size'],chunk):tasks.append((f,start,min(start+chunk,f['size'])-1))
 print('START',len(tasks),'chunks',sum(f['size'] for f in files),'bytes',flush=True)
 def fetch(task):
  f,start,end=task;error=None
  for attempt in range(1,4):
   try:
    req=urllib.request.Request(f['download_url'],headers={'Range':f'bytes={start}-{end}','User-Agent':'CAF-data-acquisition/1.0'})
    with opener.open(req,timeout=120) as resp:
     assert resp.status==206,resp.status
     assert resp.headers.get('Content-Range')==f'bytes {start}-{end}/{f["size"]}',resp.headers.get('Content-Range')
     pos=start;h=hashlib.sha256()
     while True:
      b=resp.read(1024*1024)
      if not b:break
      assert pos+len(b)<=end+1;os.pwrite(handles[f['name']],b,pos);pos+=len(b);h.update(b)
     assert pos==end+1,(pos,end+1)
    record=dict(file=f['name'],start=start,end=end,sha256=h.hexdigest(),attempt=attempt,completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
    with lock:
     with progress.open('a') as log:log.write(json.dumps(record)+'\n');log.flush()
     print('CHUNK_COMPLETE',f['name'],start,end,flush=True)
    return record
   except Exception as e:
    error=type(e).__name__+': '+str(e)[:250]
    with lock:print('CHUNK_RETRY',f['name'],start,attempt,error,flush=True)
    if attempt<3:time.sleep(3*attempt)
  raise RuntimeError(f'{f["name"]} range {start}-{end}: {error}')
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  for r in pool.map(fetch,tasks):pass
 for fd in handles.values():os.fsync(fd);os.close(fd)
 handles={};records=[]
 for f in files:
  part=D/(f['name']+'.part');md5,sha=digest_file(part);assert md5==f['computed_md5'],(f['name'],md5)
  final=D/f['name'];part.rename(final);records.append(dict(path=str(final),bytes=final.stat().st_size,md5=md5,sha256=sha,official_md5_verified=True))
  print('FILE_VERIFIED',f['name'],flush=True)
 (O/'resources/OSTEOCAR_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'status':'complete_verified','completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':records},indent=2));(O/'logs/OsteoCAR_download.exit_code').write_text('0\n')
except BaseException as e:
 for fd in handles.values():os.close(fd)
 (O/'resources/OSTEOCAR_DOWNLOAD_FAILED.json').write_text(json.dumps({'error':str(e),'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));(O/'logs/OsteoCAR_download.exit_code').write_text('1\n');traceback.print_exc();raise
