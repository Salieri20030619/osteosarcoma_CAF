
from pathlib import Path
import subprocess,os,json,datetime,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
O=Path(__file__).resolve().parents[1];R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
lock=threading.Lock();active={};done=[]
def status():
 with lock:
  f=O/'STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps({'at':datetime.datetime.now().astimezone().isoformat(),'active':dict(active),'outcomes':list(done)},indent=2));t.replace(f)
def run(name,script):
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
 with (O/'logs'/f'{name}.log').open('x') as log:
  p=subprocess.Popen([R,'-u',str(O/'code'/script)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  with lock:active[name]=p.pid
  status();rc=p.wait()
 (O/f'{name}.exit_code').write_text(str(rc))
 with lock:active.pop(name,None)
 return {'job':name,'exit_code':rc}
with ThreadPoolExecutor(max_workers=2) as pool:
 fs=[pool.submit(run,'patient_validation','patient_validation_frozen.py'),pool.submit(run,'IC_download','download_IC_frozen.py')]
 for f in as_completed(fs):
  with lock:done.append(f.result())
  status()
(O/'INITIAL_JOBS_FINISHED').write_text(datetime.datetime.now().astimezone().isoformat())
