"""Direct public SRA download with bounded retries, resumable ranges and source MD5 verification."""
from pathlib import Path
import os,sys,json,time,math,hashlib,threading,concurrent.futures,requests,fcntl,subprocess,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927')
D=O/'downloads/GSE237070_SRA';D.mkdir(parents=True,exist_ok=True)
lock=(O/'download.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
WORKERS=128;CHUNK=8*1024*1024
state_lock=threading.Lock();progress={};finished=threading.Event()
BIN=Path('/mnt/BioDisk/wangmingyuan/miniconda3/pkgs/sra-tools-3.2.1-h4304569_1/bin')
env={k:v for k,v in os.environ.items() if k.lower() not in ['http_proxy','https_proxy','all_proxy']}
env['LD_LIBRARY_PATH']='/mnt/BioDisk/wangmingyuan/miniconda3/lib'
cfg=O/'tools/sra_config';cfg.mkdir(exist_ok=True);settings=cfg/'user-settings.mkfg'
settings.write_text('/http/proxy/enabled = "false"\n');env['NCBI_SETTINGS']=str(settings)
def atomic(p,obj):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(p)
def snapshot():
 with state_lock:
  atomic(O/'DOWNLOAD_PROGRESS.json',{'time_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'workers':WORKERS,'chunk_bytes':CHUNK,'files':progress})
def monitor():
 while not finished.wait(900):snapshot()
threading.Thread(target=monitor,daemon=True).start()
def md5(path):
 h=hashlib.md5()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
 return h.hexdigest()
def process(x):
 acc=x['run'];dest=D/(acc+'.sra');part=D/(acc+'.sra.partial');journal=D/(acc+'.completed_ranges.jsonl');n=math.ceil(x['size']/CHUNK)
 with state_lock:progress[acc]={'state':'starting','size':x['size'],'gsm':x['gsm'],'patient':x['patient'],'tissue':x['tissue'],'started':time.time(),'received_this_session':0}
 snapshot()
 if not dest.exists():
  done=set()
  if journal.exists():
   for line in journal.read_text().splitlines():
    try:record=json.loads(line);assert record['source_md5']==x['md5'];done.add(int(record['chunk']))
    except (ValueError,KeyError):pass
  if part.exists():assert part.stat().st_size==x['size']
  else:assert not done;part.touch();os.truncate(part,x['size'])
  fd=os.open(part,os.O_RDWR);j=journal.open('a');jl=threading.Lock()
  with state_lock:progress[acc].update(state='downloading',completed_ranges=len(done),total_ranges=n,verified_range_bytes=sum(min(CHUNK,x['size']-i*CHUNK) for i in done))
  snapshot()
  def chunk(i):
   start=i*CHUNK;end=min(x['size'],start+CHUNK)-1;errors=[]
   for attempt in range(4):
    pos=start;t=time.time();s=requests.Session();s.trust_env=False
    try:
     with s.get(x['url'],headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(15,60)) as r:
      assert r.status_code==206,(r.status_code,r.text[:120]);assert r.headers.get('Content-Range')==f"bytes {start}-{end}/{x['size']}",dict(r.headers)
      for b in r.iter_content(65536):
       if not b:continue
       assert pos+len(b)<=end+1
       written=os.pwrite(fd,b,pos);assert written==len(b);pos+=written
       with state_lock:progress[acc]['received_this_session']+=written
       if time.time()-t>1800:raise TimeoutError('range transfer exceeded 30 minutes')
      assert pos==end+1,(pos,end)
     with jl:
      os.fsync(fd);j.write(json.dumps({'chunk':i,'start':start,'end':end,'source_md5':x['md5']})+'\n');j.flush();os.fsync(j.fileno())
     with state_lock:
      progress[acc]['completed_ranges']+=1;progress[acc]['verified_range_bytes']+=end-start+1
     return
    except Exception as e:
     errors.append(repr(e));time.sleep(min(30,2**attempt))
    finally:s.close()
   raise RuntimeError(f'{acc} chunk {i} failed: {errors}')
  failures=[]
  try:
   with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
    jobs={pool.submit(chunk,i):i for i in range(n) if i not in done}
    for future in concurrent.futures.as_completed(jobs):
     try:future.result()
     except Exception as e:failures.append(str(e))
  finally:j.close();os.close(fd)
  if failures:raise RuntimeError('; '.join(failures[:20]))
  with state_lock:progress[acc]['state']='checking_source_md5'
  snapshot();digest=md5(part);assert digest==x['md5'],f'{acc} MD5 mismatch: {digest}';part.replace(dest)
 else:
  digest=md5(dest);assert digest==x['md5'] and dest.stat().st_size==x['size']
 with state_lock:progress[acc].update(state='source_md5_passed',actual_md5=digest,bytes=dest.stat().st_size)
 snapshot()
 q=subprocess.run([str(BIN/'vdb-dump'),str(dest),'-R','1-3','-C','READ_LEN,READ_TYPE,READ_START','-f','tab'],capture_output=True,text=True,env=env,timeout=300)
 rec={'run':acc,'source_url':x['url'],'source_md5':x['md5'],'actual_md5':digest,'bytes':dest.stat().st_size,'read_probe_exit':q.returncode,'read_probe_stdout':q.stdout,'read_probe_stderr':q.stderr}
 assert q.returncode==0,rec
 assert len(q.stdout.strip().splitlines())==3,rec
 assert all('8, 28, 151' in l or '8,28,151' in l for l in q.stdout.strip().splitlines()),rec
 atomic(D/(acc+'.verified.json'),rec)
 with state_lock:progress[acc].update(state='verified',completed=time.time())
 snapshot();print(json.dumps({'event':'run_verified','run':acc,'bytes':x['size']}),flush=True)
manifest=json.loads((O/'tables/GSE237070_sra_download_manifest.json').read_text())
# A small technical run is first to validate the complete transfer/read route; all four patient-9 runs remain one biological sample.
pilot=min((x for x in manifest if x['patient']=='9'),key=lambda x:x['size'])
manifest=[pilot]+[x for x in manifest if x['run']!=pilot['run']]
atomic(O/'DOWNLOAD_STARTED.json',{'pid':os.getpid(),'time':time.time(),'total_bytes':sum(x['size'] for x in manifest),'runs':[x['run'] for x in manifest],'biological_samples':6,'patients':4,'no_proxy':True,'purpose':'source acquisition and file integrity; expression/CAF validation remains pending until matrices are reconstructed'})
code=0
try:
 for x in manifest:process(x)
 atomic(O/'DOWNLOAD_FINISHED.json',{'time':time.time(),'runs':len(manifest),'all_source_md5_and_read_layout_checks_passed':True})
except Exception as e:
 code=1;atomic(O/'DOWNLOAD_ERROR.json',{'time':time.time(),'error':repr(e),'traceback':traceback.format_exc()});print(traceback.format_exc(),flush=True)
finally:
 finished.set();snapshot();(O/'download.exit_code').write_text(str(code)+'\n')
sys.exit(code)
