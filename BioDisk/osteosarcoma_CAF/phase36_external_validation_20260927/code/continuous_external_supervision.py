
from pathlib import Path
import os,select,json,datetime,subprocess,fcntl,traceback
O=Path(__file__).resolve().parents[1];lock=(O/'external_supervisor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now().astimezone().isoformat()
def read(f):
 try:return json.loads(f.read_text())
 except Exception:return None
def state(status,**kw):
 q={'at':now(),'status':status,'monitor_pid':os.getpid(),'interval_seconds':900,'native_completion_wait':True,'initial_jobs':read(O/'STATUS.json'),'IC_download':read(O/'IC_DOWNLOAD_STATUS.json'),**kw}
 f=O/'SUPERVISION_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps(q,indent=2));t.replace(f)
 with (O/'logs/supervision_history.jsonl').open('a') as h:h.write(json.dumps(q)+'\n')
try:
 state('WAITING_FOR_IC_DOWNLOAD')
 try:
  fd=os.pidfd_open(int((O/'controller.pid').read_text()))
  while not select.select([fd],[],[],900)[0]:state('WAITING_FOR_IC_DOWNLOAD')
  os.close(fd)
 except ProcessLookupError:pass
 assert (O/'INITIAL_JOBS_FINISHED').exists() and (O/'IC_download.exit_code').read_text().strip()=='0','IC download did not complete cleanly'
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
 with (O/'logs/spatial_component_validation.log').open('x') as log:
  p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/spatial_component_validation_frozen.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  state('SPATIAL_COMPONENT_VALIDATION_RUNNING',native_pid=p.pid);rc=p.wait()
 (O/'spatial_component_validation.exit_code').write_text(str(rc));assert rc==0
 q=read(O/'delivery/SPATIAL_COMPONENT_VALIDATION_QA.json');assert q and q['samples_with_alignment']==26 and q['patient_association_tests']==32 and q['specificity_contrasts']==24
 assert (O/'SPATIAL_COMPONENT_VALIDATION_FINISHED').exists()
 state('VALIDATION_CALCULATIONS_COMPLETE_REVIEW_REQUIRED',acceptance=q)
 (O/'SUPERVISION_FINISHED').write_text(now())
except Exception:
 state('REVIEW_REQUIRED_FAILURE',error=traceback.format_exc());raise
