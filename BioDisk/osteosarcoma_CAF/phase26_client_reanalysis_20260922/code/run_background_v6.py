from common import *
import threading,fcntl,traceback,signal,select,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
H=json.loads((O/'RESERVE50_HANDOVER.json').read_text())
legacy=Path(H['backup_dir'])/'run_background.py'
exec(compile(legacy.read_text().split("\nlock=")[0],str(legacy),'exec'))
original_native_job=native_job
from memory_gate import MemoryGate,estimate_gib,available_gib,tree_pss
guard=threading.RLock()
adopted={(x['method'],x['sample']):x for x in H['adopted']}
pending_adoptions=set(adopted)
replacement_lock=(O/'runner_v6.lock').open('a')
fcntl.flock(replacement_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
legacy_lock=None

def open_pidfd(pid):
 if hasattr(os,'pidfd_open'):return os.pidfd_open(pid)
 import ctypes
 libc=ctypes.CDLL(None,use_errno=True)
 fd=libc.syscall(434,int(pid),0)
 if fd<0:raise OSError(ctypes.get_errno(),'pidfd_open')
 return fd

def proc_fields(pid):
 return (Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()

def release_legacy():
 global legacy_lock
 with guard:
  if pending_adoptions or legacy_lock is not None:return
  for c in H['controllers']:
   pid=c['pid']
   p=Path('/proc')/str(pid)
   if p.exists():
    f=proc_fields(pid)
    assert f[19]==c['start_ticks'] and f[0]=='T','Legacy controller identity/state changed'
    os.kill(pid,signal.SIGKILL) # Only controllers; all inherited native children have exited.
  legacy_lock=(O/'runner.lock').open('a')
  fcntl.flock(legacy_lock,fcntl.LOCK_EX)
  write_json(O/'LEGACY_CONTROLLERS_RETIRED.reserve50.json',{'at':now(),'controllers':H['controllers'],'all_adopted_exit_statuses_recorded':True})

def available_after_retry_reservation():
 f=O/'retry_once_20260924/MEMORY_GUARD_STATE.json'
 x=json.loads(f.read_text()) if f.exists() else {}
 extra=sum(max(0,v['budget_GiB']-(tree_pss(v['pid']) if v.get('pid') else 0)) for v in x.get('allocations',{}).values())
 return available_gib()-extra
memory_gate=MemoryGate(O/'MEMORY_GUARD_STATE.json',reserve=50,available_fn=available_after_retry_reservation)
for job in H['adopted']:
 qa=json.loads((O/'inputs'/job['sample']/'input_QA.json').read_text())
 memory_gate.tokens[job['method']+':'+job['sample']]={'budget_GiB':estimate_gib(qa),'pid':job['pid']}
with memory_gate.cond:memory_gate.status()
def launch_new_native(method,sm):
 d=O/'cnv'/method/sm;d.mkdir(parents=True,exist_ok=True)
 if (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists():return original_native_job(method,sm)
 qa=json.loads((O/'inputs'/sm/'input_QA.json').read_text())
 key=method+':'+sm
 budget=estimate_gib(qa)
 memory_gate.admit(key,budget)
 try:
  log=O/'logs'/f'{method}_{sm}.log'
  with log.open('a') as out:
   env=os.environ.copy()
   # One BLAS thread per process avoids nested BLAS x 8 worker oversubscription.
   for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='1'
   p=subprocess.Popen(['prlimit','--as=103079215104',R,'--vanilla',str(O/'code/native_cnv_cpu8.R'),method,sm],stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,env=env)
   memory_gate.attach(key,p.pid)
   rc=p.wait()
  (d/'exit_code').write_text(str(rc)+'\n')
  ok=rc==0 and (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists() and (d/'cell_results.tsv').stat().st_size>0
  return {'method':method,'status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(d),'failure':(d/'FAILED').read_text().strip() if (d/'FAILED').exists() else '', 'native_threads':8,'blas_threads_per_process':1,'admission_memory_budget_GiB':budget}
 finally:memory_gate.release(key)


def native_job(method,sm):
 key=(method,sm)
 if key not in adopted:return launch_new_native(method,sm)
 job=adopted[key];pid=job['pid'];fields=proc_fields(pid)
 assert fields[19]==job['start_ticks'] and int(fields[1])==job['parent_pid']
 fd=open_pidfd(pid)
 try:select.select([fd],[],[])
 finally:os.close(fd)
 fields=proc_fields(pid)
 assert fields[0]=='Z' and fields[19]==job['start_ticks'],'Missing inherited exit evidence'
 raw=int(fields[49]);rc=os.waitstatus_to_exitcode(raw)
 d=O/'cnv'/method/sm
 (d/'exit_code').write_text(str(rc)+'\n')
 ok=rc==0 and (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists() and (d/'cell_results.tsv').stat().st_size>0
 result={'method':method,'status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(d),'failure':(d/'FAILED').read_text().strip() if (d/'FAILED').exists() else '', 'adopted_pid':pid,'exit_status_source':'original_child_zombie_proc_stat','raw_wait_status':raw}
 write_json(d/'ADOPTED_EXIT.json',result)
 with guard:
  pending_adoptions.remove(key)
  memory_gate.release(method+':'+sm)
  release_legacy()
 return result

manifest=read(O/'tables/analysis_source_manifest.tsv').fillna('')
rows=manifest.to_dict('records')
assert len(rows)==107 and manifest['sample'].is_unique
results={r['sample']:r for r in json.loads((O/'tables/job_outcomes.json').read_text())}
assert len(results)==H['completed_samples']
groups=[g.to_dict('records') for _,g in manifest.groupby('dataset',sort=True)]
ordered=[]
while any(groups):
 for group in groups:
  if group:ordered.append(group.pop(0))
active_samples={x['sample'] for x in H['adopted']}
assert not active_samples.intersection(results)
ordered=sorted([r for r in ordered if r['sample'] not in results],key=lambda r:0 if r['sample'] in active_samples else 1)
qas={r['sample']:read_json_local(O/'inputs'/r['sample']/'input_QA.json') for r in ordered}
assert all(qas.values())
state=json.loads((O/'STATUS.json').read_text())
state.update(reserve_available_memory_GiB=50,concurrency=10,runner_pid=os.getpid(),resumed_at=now(),status='WAITING_FOR_NATIVE_CNV')
active=set()
def snapshot():
 state.update(updated_at=now(),finished=len(results),active_samples=sorted(active),queued=len(rows)-len(results)-len(active))
 write_json(O/'STATUS.json',state)
 write_json(O/'tables/job_outcomes.json',list(results.values()))
def work(row):
 sm=row['sample']
 with guard:active.add(sm);snapshot()
 try:return run_sample(row,qas[sm])
 finally:
  with guard:active.discard(sm)
def callback(f):
 sm=f.sample
 try:r=f.result()
 except Exception as e:r={'sample':sm,'status':'PROCESSING_FAILED','error':str(e),'finished':now()}
 with guard:results[sm]=r;snapshot()

if '--preflight' in sys.argv:
 import tempfile
 child=subprocess.Popen([sys.executable,'-c','raise SystemExit(7)'])
 fd=open_pidfd(child.pid);select.select([fd],[],[]);os.close(fd)
 f=proc_fields(child.pid)
 assert f[0]=='Z' and os.waitstatus_to_exitcode(int(f[49]))==7
 assert child.wait()==7
 assert len(ordered)+len(results)==107
 print(json.dumps({'preflight':'PASS','completed':len(results),'remaining':len(ordered),'adopted':len(adopted),'exit_status_capture_test':7}))
 sys.exit(0)

(O/'runner.pid').write_text(str(os.getpid())+'\n')
(O/'RUN_STARTED').write_text(str(os.getpid())+'\n')
snapshot()
write_json(O/'RESERVE50_STARTED.json',{'at':now(),'pid':os.getpid(),'concurrency':10,'preserved_completed':len(results),'adopted':H['adopted']})
release_legacy()
rc=1
try:
 with ThreadPoolExecutor(max_workers=10) as pool:
  futures=[]
  for row in ordered:
   f=pool.submit(work,row);f.sample=row['sample'];f.add_done_callback(callback);futures.append(f)
  for f in as_completed(futures):
   try:f.result()
   except Exception:pass
 assert len(results)==107
 summaries=[]
 for row in rows:
  r=results[row['sample']]
  summaries.append({'sample':row['sample'],'dataset':row['dataset'],'status':r['status'],'input_cells':r.get('input_cells'),**r.get('evidence_counts',{})})
 pd.DataFrame(summaries).to_csv(O/'tables/all107_sample_outcomes.tsv',sep='\t',index=False)
 bad=[r for r in results.values() if r['status'] in ['INPUT_FAILED','PROCESSING_FAILED']]
 qa={'finished_at':now(),'scheduled_sample_records':107,'recorded_outcomes':107,'input_or_processing_failures':len(bad),'software_failure_samples':sum(r['status']=='CNV_EVIDENCE_WITH_SOFTWARE_FAILURE' for r in results.values()),'formal_CAF_clustering_started':False,'next_action':'Review CNV and lineage identity and unresolved technical failures before CAF/myCAF clustering','all_software_success_implied':False,'concurrency':10}
 write_json(O/'delivery/CNV_STAGE_QA.json',qa)
 (O/'CNV_STAGE_FINISHED').write_text('ALL_107_SOURCE_RECORDS_HAVE_OUTCOMES\n')
 state['status']='INPUT_REPAIR_REQUIRED' if bad else 'CNV_IDENTITY_REVIEW_REQUIRED'
 rc=1 if bad else 0
 snapshot()
except Exception as e:
 state.update(status='RUNNER_FAILED',error=str(e));snapshot();traceback.print_exc()
finally:
 (O/'exit_code').write_text(str(rc)+'\n')
 (O/'RUN_FINISHED').write_text(now()+'\n')
