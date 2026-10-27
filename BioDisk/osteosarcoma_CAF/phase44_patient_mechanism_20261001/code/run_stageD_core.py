# SYNTAX_CHECK: python -m py_compile run_stageD_core.py
# RUN_COMMAND: batch_bench/bin/python -u run_stageD_core.py
# EXPECTED_SUCCESS_MARKER: STAGE_D_CORE_ACCEPTED; CORE_COMPLETE=SUCCESS; core_runner.exit_code=0
from pathlib import Path
import os,json,datetime,hashlib,subprocess,traceback
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase44_patient_mechanism_20261001';PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
state={'pid':os.getpid(),'started':datetime.datetime.now().astimezone().isoformat(),'status':'starting','sequence':['robustness_tf','clinical','prepare_QA','independent_R_QA'],'outcomes':[],'proxy_used':False,'scope':'patient stability, TF signed-response proxy and clinical proportions; final C2 candidate synthesis and communications not claimed complete'}
def write():
 f=O/'CORE_PIPELINE_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2));t.replace(f)
rc=1
try:
 for sub in ['logs','resources','delivery']:(O/sub).mkdir(exist_ok=True)
 (O/'resources/CORE_CODE_SHA256.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (O/'code').iterdir() if f.is_file()},indent=2))
 jobs=[('robustness_tf',[PY,'-u',str(O/'code/stageD_robustness_tf.py')],O/'ROBUSTNESS_TF_COMPLETE'),('clinical',[PY,'-u',str(O/'code/stageD_clinical.py')],O/'CLINICAL_COMPLETE'),('prepare_QA',[PY,'-u',str(O/'code/prepare_stageD_QA.py')],O/'CLINICAL_COMPLETE'),('independent_R_QA',[R,'--vanilla',str(O/'code/stageD_independent_QA.R')],O/'INDEPENDENT_R_QA_COMPLETE')]
 env=os.environ.copy()
 for k in list(env):
  if k.lower().endswith('_proxy'):env.pop(k,None)
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
 env['NO_PROXY']='*'
 for name,argv,marker in jobs:
  available=int(next(z for z in Path('/proc/meminfo').read_text().splitlines() if z.startswith('MemAvailable:')).split()[1])/1048576
  if available<50:raise RuntimeError('Memory reserve below 50 GiB; job not launched')
  with (O/'logs'/f'{name}.log').open('x') as log:
   child=subprocess.Popen(['/usr/bin/time','-v',*argv],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env);state.update(status='running',active_job=name,active_pid=child.pid);write();print('STARTED',name,child.pid,flush=True);result=child.wait()
  (O/f'{name}.runner_exit_code').write_text(str(result)+'\n');state['outcomes'].append({'job':name,'exit_code':result,'marker':str(marker),'marker_present':marker.exists()});write()
  if result!=0 or not marker.exists() or marker.read_text().strip()!='SUCCESS':raise RuntimeError(f'{name} failed acceptance')
  print('ACCEPTED',name,flush=True)
 state.update(status='complete',completed=datetime.datetime.now().astimezone().isoformat(),active_job=None,active_pid=None);write();(O/'CORE_COMPLETE').write_text('SUCCESS\n');print('STAGE_D_CORE_ACCEPTED',flush=True);rc=0
except Exception:traceback.print_exc();state.update(status='failed',error=traceback.format_exc());write();(O/'CORE_FAILED').write_text(traceback.format_exc())
finally:(O/'core_runner.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)