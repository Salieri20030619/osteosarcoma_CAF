# SYNTAX_CHECK: python -m py_compile stageC_controller.py
# RUN_COMMAND: python start_stageC.py
# EXPECTED_SUCCESS_MARKER: SEQUENTIAL_STAGE_C_COMPLETE, PIPELINE_COMPLETE=SUCCESS, controller.exit_code=0
from pathlib import Path
import os,json,datetime,subprocess,hashlib,traceback,resource
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase43_full_transcriptome_discovery_20261001';PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
for sub in ['logs','code','delivery','resources']:(O/sub).mkdir(parents=True,exist_ok=True)
state={'pid':os.getpid(),'started':datetime.datetime.now().astimezone().isoformat(),'status':'starting','sequence':['full_gene_statistics','prepare_independent_QA','independent_R_QA','open_pathways'],'outcomes':[],'proxy_used':False,'memory_reserve_GiB':50,'threads':4}
def status():
 f=O/'PIPELINE_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2));t.replace(f)
def run(name,argv,marker):
 avail=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])/1048576
 if avail<50:raise RuntimeError(f'Available memory {avail:.1f} GiB below preserved 50 GiB; no job launched')
 env=os.environ.copy()
 for key in list(env):
  if key.lower().endswith('_proxy'):env.pop(key,None)
 for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[key]='4'
 env['NO_PROXY']='*'
 with (O/'logs'/f'{name}.log').open('x') as log:
  child=subprocess.Popen(['/usr/bin/time','-v',*argv],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  state.update(status='running',active_job=name,active_pid=child.pid);(O/f'{name}.pid').write_text(str(child.pid)+'\n');status();rc=child.wait()
 (O/f'{name}.controller_exit_code').write_text(str(rc)+'\n');state['outcomes'].append({'job':name,'exit_code':rc,'expected_marker':str(marker),'marker_present':marker.exists()});status()
 if rc!=0 or not marker.exists() or marker.read_text().strip()!='SUCCESS':raise RuntimeError(f'{name} failed acceptance, exit={rc}, marker={marker}')
rc=1
try:
 for phase,exitname in [('phase41_full_transcriptome_setup_20261001','setup.exit_code'),('phase42_F3F4_identity_states_20261001','identity_states.exit_code')]:
  q=P/phase;assert (q/exitname).read_text().strip()=='0' and (q/'COMPLETE').read_text().strip()=='SUCCESS'
  failed=q/'FAILED'
  if failed.exists():failed.rename(q/'logs/FAILED_historical_before_success.txt')
 run('full_gene_statistics',[PY,'-u',str(O/'code/stageC_full_gene_statistics.py')],O/'C1_COMPLETE')
 run('prepare_independent_QA',[PY,'-u',str(O/'code/prepare_independent_QA.py')],O/'C1_COMPLETE')
 run('independent_R_QA',[R,'--vanilla',str(O/'code/stageC_independent_QA.R')],O/'INDEPENDENT_R_QA_COMPLETE')
 run('open_pathways',[PY,'-u',str(O/'code/stageC_pathways.py')],O/'C2_COMPLETE')
 state.update(status='complete',completed=datetime.datetime.now().astimezone().isoformat(),active_job=None,active_pid=None);status();(O/'PIPELINE_COMPLETE').write_text('SUCCESS\n');print('SEQUENTIAL_STAGE_C_COMPLETE',flush=True);rc=0
except Exception:
 traceback.print_exc();state.update(status='failed',error=traceback.format_exc());status();(O/'PIPELINE_FAILED').write_text(traceback.format_exc())
finally:(O/'controller.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)