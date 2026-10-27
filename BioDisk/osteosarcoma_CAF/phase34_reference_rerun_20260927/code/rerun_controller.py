
from pathlib import Path
import os,json,subprocess,datetime,threading,hashlib,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed
O=Path(__file__).resolve().parents[1]
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
protocol=json.loads((O/'resources/RERUN_PROTOCOL.json').read_text())
lock=threading.Lock();active={};outcomes=[]
def now():return datetime.datetime.now().astimezone().isoformat()
def state():
 with lock:
  f=O/'STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps({'at':now(),'status':'RUNNING' if len(outcomes)<4 else 'COMPLETE' if all(x['accepted'] for x in outcomes) else 'FAILED','scheduled':4,'finished':len(outcomes),'active':dict(active),'outcomes':list(outcomes)},indent=2));t.replace(f)
def run(j):
 sm=j['sample'];v=j['variant'];d=O/'cnv_reference'/sm/v;d.mkdir(parents=True)
 script=O/'code/reference_sensitivity_frozen.R';wrapper=O/'code/run_frozen.R'
 assert hashlib.sha256(script.read_bytes()).hexdigest()==protocol['script_sha256']
 assert hashlib.sha256(wrapper.read_bytes()).hexdigest()==protocol['wrapper_sha256']
 env=os.environ.copy();env['CAF_FROZEN_R_SCRIPT']=str(script)
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
 with (O/'logs'/f'{sm}__{v}.log').open('x') as log:
  p=subprocess.Popen([R,'--vanilla',str(wrapper),sm,v],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  with lock:active[sm+'__'+v]=p.pid
  state();rc=p.wait()
 (d/'exit_code').write_text(str(rc))
 import pandas as pd
 accepted=False;detail=''
 try:
  assert rc==0 and (d/'COMPLETE').exists() and not (d/'FAILED').exists()
  x=pd.read_csv(d/'cell_results.tsv',sep='\t');a=pd.read_csv(O/'inputs'/sm/f'{v}_annotations.tsv',sep='\t',header=None)
  assert len(x)==j['input_cells'] and x.cell_id.is_unique and len(a)==len(x)
  assert set(x.call)<= {'alteration_called','no_alteration_called','not_returned_by_software'}
  assert (d/'run.final.infercnv_obj').stat().st_size>0
  assert hashlib.sha256(script.read_bytes()).hexdigest()==protocol['script_sha256']
  accepted=True
 except Exception as e:detail=repr(e)
 with lock:active.pop(sm+'__'+v,None)
 return dict(sample=sm,variant=v,exit_code=rc,accepted=accepted,detail=detail,script_unchanged=hashlib.sha256(script.read_bytes()).hexdigest()==protocol['script_sha256'])
state()
try:
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(run,j) for j in protocol['jobs']]):
   with lock:outcomes.append(f.result())
   state()
 (O/'exit_code').write_text('0' if all(x['accepted'] for x in outcomes) else '1')
 (O/'FINISHED').write_text(now());state()
except Exception:
 (O/'exit_code').write_text('1');traceback.print_exc();raise
