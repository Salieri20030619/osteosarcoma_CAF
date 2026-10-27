
from pathlib import Path
import os,json,subprocess,datetime,threading,time,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase33_client_report_20260926')
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
jobs=[(sm,v) for sm in ['SCPCS000614','SCPCS000610'] for v in ['source_reference','lymphoid_reference']]
lock=threading.Lock();active={};outcomes=[]
def now():return datetime.datetime.now().astimezone().isoformat()
def state():
 with lock:
  f=O/'REFERENCE_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps({'at':now(),'scheduled':4,'finished':len(outcomes),'active':dict(active),'outcomes':outcomes},indent=2));t.replace(f)
def run(job):
 sm,var=job;d=O/'cnv_reference'/sm/var;d.mkdir(parents=True,exist_ok=True)
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
 with (O/'logs'/f'cnv_{sm}_{var}.log').open('x') as log:
  p=subprocess.Popen([R,'--vanilla',str(O/'code/reference_sensitivity.R'),sm,var],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  with lock:active[sm+'__'+var]=p.pid
  state();rc=p.wait()
 (d/'exit_code').write_text(str(rc))
 ok=rc==0 and (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists()
 with lock:active.pop(sm+'__'+var,None)
 return {'sample':sm,'variant':var,'exit_code':rc,'complete':ok,'failure':(d/'FAILED').read_text() if (d/'FAILED').exists() else ''}
try:
 with ThreadPoolExecutor(max_workers=2) as pool:
  fs=[pool.submit(run,j) for j in jobs]
  for f in as_completed(fs):
   with lock:outcomes.append(f.result())
   state()
 (O/'reference.exit_code').write_text('0' if all(r['complete'] for r in outcomes) else '1')
 (O/'REFERENCE_FINISHED').write_text(now());state()
except Exception:
 traceback.print_exc();(O/'reference.exit_code').write_text('1');raise
