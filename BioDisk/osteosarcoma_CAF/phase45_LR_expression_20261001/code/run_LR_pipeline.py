from pathlib import Path
import os,sys,subprocess,json,datetime,traceback,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase45_LR_expression_20261001');PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
for sub in ['logs','delivery','resources']:(O/sub).mkdir(exist_ok=True)
if '--start' in sys.argv:
 assert not (O/'controller.pid').exists(),'already launched'
 env=os.environ.copy()
 for k in list(env):
  if k.lower().endswith('_proxy'):env.pop(k,None)
 env['NO_PROXY']='*'
 with (O/'logs/controller.log').open('x') as log:p=subprocess.Popen([PY,'-u',__file__],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
 (O/'controller.pid').write_text(str(p.pid)+'\n');print('PATIENT_LR_PIPELINE_STARTED',p.pid);raise SystemExit(0)
state={'pid':os.getpid(),'status':'running','started':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'outcomes':[]};rc=1
try:
 (O/'resources/CODE_SHA256.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (O/'code').iterdir() if f.is_file()},indent=2))
 for name,argv,marker in [('LR_expression',[PY,'-u',str(O/'code/stageD_LR_expression.py')],O/'COMPLETE'),('independent_R',[R,'--vanilla',str(O/'code/stageD_LR_QA.R')],O/'INDEPENDENT_R_QA_COMPLETE')]:
  mem=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])/1048576
  if mem<50:raise RuntimeError('memory reserve below 50 GiB')
  env=os.environ.copy()
  for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
  with (O/'logs'/f'{name}.log').open('x') as log:
   p=subprocess.Popen(['/usr/bin/time','-v',*argv],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env);state.update(active_job=name,active_pid=p.pid);(O/'PIPELINE_STATUS.json').write_text(json.dumps(state,indent=2));code=p.wait()
  (O/f'{name}.runner_exit_code').write_text(str(code)+'\n');state['outcomes'].append({'job':name,'exit_code':code,'marker_present':marker.exists()});(O/'PIPELINE_STATUS.json').write_text(json.dumps(state,indent=2))
  assert code==0 and marker.read_text().strip()=='SUCCESS',name
 state.update(status='complete',completed=datetime.datetime.now().astimezone().isoformat(),active_job=None,active_pid=None);(O/'PIPELINE_STATUS.json').write_text(json.dumps(state,indent=2));(O/'PIPELINE_COMPLETE').write_text('SUCCESS\n');print('PATIENT_LR_PIPELINE_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();state.update(status='failed',error=traceback.format_exc());(O/'PIPELINE_STATUS.json').write_text(json.dumps(state,indent=2));(O/'PIPELINE_FAILED').write_text(traceback.format_exc())
finally:(O/'controller.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)