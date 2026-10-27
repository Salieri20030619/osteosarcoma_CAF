"""SYNTAX_CHECK: python -m py_compile code/run_cohort.py
RUN_COMMAND: python -u code/run_cohort.py
EXPECTED_SUCCESS_MARKER: PHASE6_COHORT_AND_INTEGRATION_FINISHED
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,os,subprocess,sys,datetime,csv,fcntl
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase6_nonpdx_cnv_20260916')
lock=(p/'logs/cohort.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert (p/'INPUTS_READY').exists()
rows=list(csv.DictReader((p/'review/input_summary.tsv').open(),delimiter='\t'))
rows=sorted([r for r in rows if r['cnv_eligible']=='True'],key=lambda r:int(r['run_cells']))
assert rows
(p/'logs/cohort.pid').write_text(str(os.getpid()))
state={'state':'DUAL_CNV_RUNNING','started_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'eligible_samples':len(rows),'planned_method_jobs':2*len(rows),'parallel_method_lanes':2,'threads_per_method_job':4,'virtual_memory_limit_per_job_GiB':64,'completion_wait_policy':'server-side process monitoring at 900-second intervals'}
(p/'STATUS.json').write_text(json.dumps(state,indent=2))
def lane(method):
 result=[]
 for row in rows:
  s=row['sample'];print('JOB_START',method,s,row['run_cells'],flush=True)
  with (p/'logs'/(method+'_'+s+'.log')).open('a') as f:
   rc=subprocess.run(['bash',str(p/'code/run_cnv_job.sh'),method,s],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT).returncode
  done=(p/method/s/'COMPLETE').exists()
  result.append({'method':method,'sample':s,'exit_code':rc,'complete_marker':done})
  (p/'logs'/(method+'_cohort_jobs.json')).write_text(json.dumps(result,indent=2))
  print('JOB_EXIT',method,s,rc,done,flush=True)
 return result
with ThreadPoolExecutor(max_workers=2) as pool:
 results=list(pool.map(lane,['infercnv','copykat']))
flat=sum(results,[]);(p/'logs/cohort_jobs.json').write_text(json.dumps(flat,indent=2))
print('ALL_METHOD_JOBS_EXITED',len(flat),flush=True)
with (p/'logs/integrate.log').open('a') as f:rc=subprocess.run([sys.executable,'-u',str(p/'code/integrate_identity.py')],stdout=f,stderr=subprocess.STDOUT).returncode
code=0 if rc==0 and all(x['exit_code']==0 and x['complete_marker'] for x in flat) else 2
(p/'logs/cohort.exit_code').write_text(str(code))
(p/'logs/cohort.finished_at').write_text(datetime.datetime.now(datetime.timezone.utc).isoformat())
print('PHASE6_COHORT_AND_INTEGRATION_FINISHED',code,flush=True)
sys.exit(code)
