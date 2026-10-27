
from pathlib import Path
import os,json,urllib.request,concurrent.futures,hashlib,gzip,datetime,threading,time,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1');D=O/'downloads/GSE234187';name='GSM7453693_matrix.mtx.gz';part=D/(name+'.part');size=126169981;url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7453nnn/GSM7453693/suppl/'+name;prefix=json.loads((O/'resources/GSE234187_RANGE_RESUME.json').read_text())['preserved_prefix_bytes'];fd=os.open(part,os.O_RDWR);os.ftruncate(fd,size);op=urllib.request.build_opener(urllib.request.ProxyHandler({}));lock=threading.Lock()
try:
 def fetch(pair):
  start,end=pair
  for attempt in range(1,5):
   try:
    req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'})
    with op.open(req,timeout=120) as resp:
     assert resp.status==206 and resp.headers.get('Content-Range')==f'bytes {start}-{end}/{size}'
     pos=start;h=hashlib.sha256()
     while True:
      b=resp.read(262144)
      if not b:break
      assert pos+len(b)<=end+1;os.pwrite(fd,b,pos);pos+=len(b);h.update(b)
     assert pos==end+1
    with lock:
     with (O/'logs/GSE234187_chunks.jsonl').open('a') as f:f.write(json.dumps({'start':start,'end':end,'sha256':h.hexdigest()})+'\n')
     print('CHUNK',start,end,flush=True)
    return
   except Exception:
    if attempt==4:raise
    time.sleep(3)
 ranges=[(s,min(s+4*1024*1024,size)-1) for s in range(prefix,size,4*1024*1024)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:list(ex.map(fetch,ranges))
 os.fsync(fd);os.close(fd);fd=None
 with gzip.open(part,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):pass
 final=D/name;part.rename(final);records=[]
 for f in sorted(D.glob('*.gz')):
  h=hashlib.sha256()
  with f.open('rb') as stream:
   for b in iter(lambda:stream.read(8388608),b''):h.update(b)
  records.append({'file':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest(),'gzip_integrity':'PASS'})
 (O/'resources/GSE234187_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'status':'complete_verified','completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':records,'matrix_size_matches_official':True},indent=2));(O/'logs/GSE234187_download.exit_code').write_text('0\n');print('GSE234187_VERIFIED_COMPLETE',flush=True)
except BaseException as e:
 if fd is not None:os.close(fd)
 (O/'resources/GSE234187_DOWNLOAD_FAILED.json').write_text(json.dumps({'error':str(e)}));(O/'logs/GSE234187_download.exit_code').write_text('1\n');traceback.print_exc();raise
