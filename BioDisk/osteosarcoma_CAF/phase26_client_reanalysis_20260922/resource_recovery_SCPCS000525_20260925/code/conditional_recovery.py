from pathlib import Path
import os,sys,json,datetime,fcntl,select,subprocess,resource,signal,hashlib,traceback
BASE=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase26_client_reanalysis_20260922')
PREVIOUS=BASE/'retry_once_20260924'
OUT=BASE/'resource_recovery_SCPCS000525_20260925'
SM='SCPCS000525'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
sys.path.insert(0,str(PREVIOUS/'code'))
import retry_early_v3 as old
sys.path.insert(0,str(BASE/'code'))
from memory_gate import MemoryGate,available_gib,tree_pss,estimate_gib

def now():return datetime.datetime.now().astimezone().isoformat()
def read(p):return json.loads(p.read_text())
def write(p,x):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decision(job,exit_code,complete,nonempty_result,failed_marker):
 if job['status'].startswith('complete') and job.get('exit_code')==0 and exit_code==0 and complete and nonempty_result:return 'KEEP_SUCCESS'
 if job.get('exit_code')==exit_code and job.get('status')=='software_failed' and (exit_code!=0 or failed_marker):return 'RERUN_FAILED_SAMPLE'
 raise RuntimeError('Inconsistent previous sample outcome; refuse automatic rerun')
def finish(state,status,qa):
 write(OUT/'delivery/QA.json',qa)
 state.update(status=status,at=now(),native_pid=None);write(OUT/'STATUS.json',state)
 (OUT/'RUN_FINISHED').write_text(now()+'\n');(OUT/'exit_code').write_text('0\n')
def main():
 lock=(OUT/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (OUT/'RUN_FINISHED').exists()
 plan=read(OUT/'PLAN.json');pid=plan['previous_controller_pid']
 state={'status':'WAITING_FOR_CURRENT_RETRY_COMPLETION','at':now(),'pid':os.getpid(),'sample':SM,'previous_controller_pid':pid,'native_pid':None}
 write(OUT/'STATUS.json',state);(OUT/'runner.pid').write_text(str(os.getpid())+'\n')
 try:
  fields=old.pid_fields(pid);assert fields[19]==plan['previous_controller_start_ticks']
  fd=old.pidfd(pid)
  try:select.select([fd],[],[])
  finally:os.close(fd)
 except FileNotFoundError:pass
 assert int((PREVIOUS/'exit_code').read_text())==0,'Previous controller did not complete cleanly'
 assert (PREVIOUS/'RETRY_STAGE_FINISHED').exists()
 previous_qa=read(PREVIOUS/'delivery/RETRY_QA.json')
 assert previous_qa['first_pass_outcomes']==107
 rows=read(PREVIOUS/'tables/job_outcomes_after_retry.json')
 assert len(rows)==len({r['sample'] for r in rows})==107
 rr=next(r for r in rows if r['sample']==SM)
 job=next(j for j in rr['jobs'] if j['method']=='copykat');d=Path(job['path'])
 rc=int((d/'exit_code').read_text())
 dec=decision(job,rc,(d/'COMPLETE').exists(),(d/'cell_results.tsv').exists() and (d/'cell_results.tsv').stat().st_size>0,(d/'FAILED').exists())
 record={'at':now(),'decision':dec,'previous_job':job,'previous_QA':previous_qa};write(OUT/'PREVIOUS_OUTCOME.json',record)
 if dec=='KEEP_SUCCESS':
  n=old.validate_cells(d/'cell_results.tsv',SM)
  finish(state,'CURRENT_RETRY_SUCCEEDED_NO_RESTART',{'at':now(),'sample':SM,'validated_cells':n,'restarted':False,'original_result_retained':True});return
 state.update(status='WAITING_FOR_MEMORY',at=now(),previous_failure=job.get('failure'));write(OUT/'STATUS.json',state)
 assert resource.getrlimit(resource.RLIMIT_AS)==(-1,-1),'Recovery controller must have unlimited address space'
 gate=MemoryGate(OUT/'MEMORY_GUARD_STATE.json',reserve=50)
 key='copykat:'+SM;budget=estimate_gib(read(BASE/'inputs'/SM/'input_QA.json'));gate.admit(key,budget)
 original_files=[BASE/'tables/job_outcomes.json',PREVIOUS/'retry_outcomes.json',PREVIOUS/'tables/job_outcomes_after_retry.json']
 hashes={str(f):sha(f) for f in original_files}
 native=OUT/'code/native_cnv_recovery.R';assert sha(native)==plan['native_script_sha256']
 result_dir=OUT/'cnv/copykat'/SM;result_dir.mkdir(parents=True,exist_ok=False)
 env=os.environ.copy();env['CAF_RETRY_ROOT']=str(OUT)
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='1'
 stopped_reason='';validated=None
 (result_dir/'ATTEMPT_STARTED').write_text(now()+'\n')
 try:
  with (OUT/'logs/copykat_SCPCS000525.log').open('x') as log:
   child=subprocess.Popen([R,'--vanilla',str(native),'copykat',SM],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
   gate.attach(key,child.pid)
   state.update(status='RUNNING_UNCAPPED_RECOVERY',native_pid=child.pid,at=now(),address_space_limit_GiB=None,reserve_GiB=50);write(OUT/'STATUS.json',state)
   assert resource.prlimit(child.pid,resource.RLIMIT_AS)==(-1,-1)
   write(OUT/'NATIVE_START.json',{'at':now(),'pid':child.pid,'start_ticks':old.pid_fields(child.pid)[19],'address_space_limits':list(resource.prlimit(child.pid,resource.RLIMIT_AS)),'memory_budget_GiB':budget})
   fd=old.pidfd(child.pid)
   try:
    while not select.select([fd],[],[],30)[0]:
     available=available_gib()
     write(OUT/'ACTIVE_MEMORY_GUARD.json',{'at':now(),'available_GiB':available,'minimum_available_GiB':50,'pid':child.pid,'action':'stop' if available<50 else 'continue'})
     if available<50:
      stopped_reason='Stopped by resource guard: host available memory below 50 GiB'
      os.killpg(child.pid,signal.SIGTERM)
      if not select.select([fd],[],[],15)[0]:os.killpg(child.pid,signal.SIGKILL)
      break
     with gate.cond:gate.status()
     state['at']=now();write(OUT/'STATUS.json',state)
   finally:os.close(fd)
   rc=child.wait()
  (result_dir/'exit_code').write_text(str(rc)+'\n')
  ok=rc==0 and (result_dir/'COMPLETE').exists() and (result_dir/'cell_results.tsv').exists()
  validation_error=''
  if ok:
   try:validated=old.validate_cells(result_dir/'cell_results.tsv',SM)
   except Exception as e:ok=False;validation_error=str(e)
 finally:gate.release(key)
 failure=stopped_reason or ((result_dir/'FAILED').read_text().strip() if (result_dir/'FAILED').exists() else validation_error)
 new={'method':'copykat','status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(result_dir),'failure':failure,'validated_cells':validated,'native_threads':8,'address_space_limit_GiB':None,'resource_recovery_attempt':1,'previous_capped_retry':str(d)}
 rr['jobs']=[new if j['method']=='copykat' else j for j in rr['jobs']]
 rr['status']='CNV_EVIDENCE_READY' if all(old.success(j) for j in rr['jobs']) else 'CNV_EVIDENCE_WITH_SOFTWARE_FAILURE'
 ns={};src=(BASE/'code/run_background.py').read_text().split('\nlock=')[0];exec(compile(src,str(BASE/'code/run_background.py'),'exec'),ns);ns['O']=OUT
 manifest=ns['read'](BASE/'tables/analysis_source_manifest.tsv').fillna('');row=next(x for x in manifest.to_dict('records') if x['sample']==SM)
 rr['evidence_counts']=ns['combine'](row,rr['jobs']);rr['resource_recovery_finished_at']=now()
 write(OUT/'tables/all107_outcomes_after_resource_recovery.json',rows)
 assert all(sha(Path(f))==v for f,v in hashes.items()),'Previous outputs changed'
 before={r['sample']:r for r in read(PREVIOUS/'tables/job_outcomes_after_retry.json')}
 assert all(r==before[r['sample']] for r in rows if r['sample']!=SM)
 qa={'at':now(),'sample':SM,'restarted':True,'native_exit_code':rc,'native_success':ok,'validated_cells':validated,'failure':failure,'merged_sample_records':107,'unchanged_other_samples':106,'previous_outputs_unchanged':True,'formal_CAF_clustering_started':False,'additional_automatic_restarts':False}
 finish(state,'RESOURCE_RECOVERY_COMPLETE_REVIEW_REQUIRED' if ok else 'RESOURCE_RECOVERY_FAILED_REVIEW_REQUIRED',qa)
if __name__=='__main__':
 try:main()
 except Exception as e:
  traceback.print_exc();write(OUT/'ERROR.json',{'at':now(),'error':str(e)});write(OUT/'STATUS.json',{'status':'CONTROLLER_FAILED','at':now(),'error':str(e)});(OUT/'exit_code').write_text('1\n');sys.exit(1)
