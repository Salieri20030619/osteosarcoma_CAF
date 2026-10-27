from pathlib import Path
import os,json,datetime,subprocess,fcntl,traceback,threading,time
from concurrent.futures import ThreadPoolExecutor,as_completed
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase31_customer_completion_20260926';B=P/'phase26_client_reanalysis_20260922'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):
 t=f.with_name(f.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(f)
guard=threading.Lock();active={};results=[];admit=threading.Condition()
def available():
 return int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))/1024**2
def status(s):
 with guard:writej(O/'METHODS_STATUS.json',{'status':s,'at':now(),'scheduled':len(jobs),'finished':len(results),'active':dict(active),'max_R_jobs':2,'threads_per_R_job':4,'available_memory_GiB':available(),'reserve_available_GiB':50})
def run(job):
 kind,sm=job;d=O/kind/sm;d.mkdir(parents=True,exist_ok=True)
 if (d/'COMPLETE').exists() and (d/'exit_code').read_text().strip()=='0':return {'sample':sm,'method':kind,'status':'complete_cached','exit_code':0,'path':str(d)}
 # Resource headroom is not a biological sample gate; wait for admission rather than dropping samples.
 with admit:
  while available()<50+32:
   with guard:active[kind+':'+sm]={'state':'WAITING_FOR_MEMORY'}
   status('RUNNING');admit.wait(timeout=60)
  with guard:active[kind+':'+sm]={'state':'STARTING'}
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
 with (O/'logs'/f'{kind}_{sm}.log').open('x') as log:
  p=subprocess.Popen([R,'--vanilla',str(O/'code/native_quality.R'),kind,sm],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  with guard:active[kind+':'+sm]={'state':'RUNNING','pid':p.pid}
  status('RUNNING');rc=p.wait()
 (d/'exit_code').write_text(str(rc)+chr(10));ok=rc==0 and (d/'COMPLETE').exists()
 if ok:
  expected='corrected_counts.mtx.gz' if kind=='decontx' else 'doublet_results.tsv'
  ok=(d/expected).exists() and (d/expected).stat().st_size>0
 result={'sample':sm,'method':kind,'status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(d),'failure':(d/'FAILED').read_text() if (d/'FAILED').exists() else ''}
 with guard:active.pop(kind+':'+sm,None)
 with admit:admit.notify_all()
 return result
lock=(O/'methods.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
plan=json.loads((O/'METHOD_JOBS.json').read_text());jobs=[(r['method'],r['sample']) for r in plan]
rc=1
try:
 status('RUNNING')
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures={pool.submit(run,j):j for j in jobs}
  for f in as_completed(futures):
   try:r=f.result()
   except Exception as ex:
    kind,sm=futures[f];r={'sample':sm,'method':kind,'status':'controller_job_error','error':str(ex)}
   with guard:
    results.append(r);writej(O/'tables/method_outcomes.json',results)
   status('RUNNING');print('FINISHED',r,flush=True)
 qa={'at':now(),'scheduled':len(jobs),'recorded':len(results),'decontx_scheduled':sum(k=='decontx' for k,s in jobs),'complete':sum(r['status'].startswith('complete') for r in results),'failed':sum(not r['status'].startswith('complete') for r in results),'no_cells_automatically_dropped':True}
 assert len(results)==len(jobs);writej(O/'delivery/METHODS_QA.json',qa)
 (O/'METHODS_FINISHED').write_text(now()+chr(10));rc=0 if qa['failed']==0 else 1;status('ALL_METHOD_JOBS_HAVE_OUTCOMES')
except Exception as ex:writej(O/'METHODS_ERROR.json',{'at':now(),'error':str(ex)});traceback.print_exc()
finally:(O/'methods.exit_code').write_text(str(rc)+chr(10))
raise SystemExit(rc)
