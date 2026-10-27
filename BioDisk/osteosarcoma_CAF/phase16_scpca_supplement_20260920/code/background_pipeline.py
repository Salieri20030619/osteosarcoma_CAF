"""Persistent server-side orchestration; native subprocess waits, no interactive dependency."""
from pathlib import Path
import os,json,subprocess,datetime,fcntl,traceback,time,csv,threading
from concurrent.futures import ThreadPoolExecutor
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase16_scpca_supplement_20260920'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
mutex=threading.Lock();records=[]
state={'status':'PREPARING_FOLLOWUP','started_utc':now(),'pid':os.getpid(),'automatic_chaining':True,'steps':['prepare','ambient_and_dual_CNV','identity_panels_paired_surface','report_integrity'],'active_jobs':{},'completed_jobs':0}
def dump(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2));t.replace(p)
def save_state():dump(O/'BACKGROUND_STATUS.json',state)
def project_state(status):
 p=P/'PROJECT_STATE.json';v=json.loads(p.read_text());v.update(current_stage=O.name,current_status=status,status=status,dynamic_status_file=str(O/'BACKGROUND_STATUS.json'),updated_at=now());v['phase16_background']=state.copy()
 if (O/'ANALYSIS_COMPLETE').exists():v['latest_report']=str(O/'delivery/补充分析结果.md');v['next_action']='Review completed ScPCA evidence and customer goal assessment; experimental validation remains external'
 dump(p,v)
def run_stage(name):
 with (O/'logs'/f'{name}.log').open('a') as f:rc=subprocess.run([PY,'-u',str(O/'code'/f'{name}.py')],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT).returncode
 (O/'logs'/f'{name}.exit_code').write_text(str(rc)+'\n');assert rc==0,(name,rc)
def mem_available_kib():
 for line in Path('/proc/meminfo').read_text().splitlines():
  if line.startswith('MemAvailable:'):return int(line.split()[1])
 raise RuntimeError('Cannot read available memory')
def validate(method,sm):
 d=O/method/sm;assert (d/'COMPLETE').exists(), 'completion marker absent'
 if method=='decontx':
  import pandas as pd
  inp=pd.read_csv(O/'inputs'/sm/'cells.tsv',sep='\t');v=pd.read_csv(d/'cell_contamination.tsv',sep='\t')
  assert v.cell_id.tolist()==inp.cell_id.tolist() and v.contamination.between(0,1).all()
  assert (d/'corrected_counts.mtx.gz').stat().st_size>0
 else:
  inp=list(csv.DictReader((O/'cnv_inputs'/sm/'cells.tsv').open(),delimiter='\t'));v=list(csv.DictReader((d/'cell_results.tsv').open(),delimiter='\t'))
  assert [(x['cell_id'],x['tool_id']) for x in v]==[(x['cell_id'],x['tool_id']) for x in inp]
  par=json.loads((d/'parameters.json').read_text());assert par['reference_cells']>=30 and par['heldout_cells']>=20
def lane(method,rows):
 for row in rows:
  sm=row['sample'];key=method+'_'+sm
  with mutex:state['active_jobs'][method]={'sample':sm,'started_utc':now(),'state':'waiting_resource'};save_state()
  # Persistent native monitor, 15-minute resource checks only when needed.
  while mem_available_kib()<100*1024*1024:time.sleep(900)
  command=[R,'--vanilla',str(O/'code/run_decontx.R'),sm] if method=='decontx' else ['bash',str(O/'code/run_cnv_job.sh'),method,sm]
  attempts=[];passed=False;why=''
  for attempt in range(2):
   with mutex:state['active_jobs'][method].update(state='running',attempt=attempt+1);save_state()
   with (O/'logs'/f'{key}.log').open('a') as f:rc=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4','MKL_NUM_THREADS':'4'}).returncode
   (O/'logs'/f'{key}.exit_code').write_text(str(rc)+'\n');attempts.append(rc)
   try:
    assert rc==0, f'exit {rc}'
    validate(method,sm);passed=True;why='';break
   except Exception as e:why=repr(e)
  rec={'method':method,'sample':sm,'passed':passed,'attempt_exit_codes':attempts,'error':why,'finished_utc':now()}
  with mutex:
   records.append(rec);dump(O/'logs/followup_jobs.json',records);state['completed_jobs']=len(records);state['active_jobs'].pop(method,None);save_state()
  print('JOB_FINISHED',json.dumps(rec),flush=True)
def main():
 lock=(O/'logs/background.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 if (O/'ANALYSIS_COMPLETE').exists():raise RuntimeError('Analysis already complete; refuse duplicate launch')
 save_state();project_state(state['status']);run_stage('prepare_followup')
 assert (O/'FOLLOWUP_INPUTS_READY').exists()
 rows=list(csv.DictReader((O/'review/input_summary.tsv').open(),delimiter='\t'))
 amb=[r for r in rows if r['ambient_eligible']=='True'];cnv=[r for r in rows if r['cnv_eligible']=='True']
 amb.sort(key=lambda r:int(r['cells']));cnv.sort(key=lambda r:int(r['cells']))
 state.update(status='AMBIENT_AND_DUAL_CNV_RUNNING',ambient_samples=len(amb),CNV_samples=len(cnv),planned_jobs=len(amb)+2*len(cnv),data_unavailable=[r for r in rows if r['cnv_eligible']!='True' or r['ambient_eligible']!='True']);save_state();project_state(state['status'])
 with ThreadPoolExecutor(max_workers=3) as pool:
  futures=[pool.submit(lane,'decontx',amb),pool.submit(lane,'infercnv',cnv),pool.submit(lane,'copykat',cnv)]
  for future in futures:future.result()
 assert len(records)==state['planned_jobs'] and all(r['passed'] for r in records),'One or more methods failed; downstream acceptance withheld'
 state['status']='IDENTITY_AND_CLIENT_ANALYSIS_RUNNING';save_state();project_state(state['status']);run_stage('finalize_followup')
 assert (O/'ANALYSIS_COMPLETE').exists()
 qa=json.loads((O/'delivery/QA.json').read_text());assert qa['all_eligible_computations_complete'] and qa['reviewed_cells']==44242
 state.update(status='ALL_ELIGIBLE_ANALYSES_COMPLETE',finished_utc=now(),QA=qa,report=str(O/'delivery/补充分析结果.md'));save_state();project_state(state['status'])
 (O/'logs/background.exit_code').write_text('0\n');print('BACKGROUND_PIPELINE_COMPLETE',flush=True)
try:main()
except Exception:
 state.update(status='BACKGROUND_PIPELINE_FAILED',failed_utc=now(),error=traceback.format_exc());save_state();project_state(state['status']);(O/'logs/background.exit_code').write_text('1\n');raise
