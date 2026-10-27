from pathlib import Path
import os,select,json,subprocess,datetime,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');BASE=P/'phase90_patient_validation_20261015'
lock=(BASE/'detached_review_events.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
configs=[('NCC_all',BASE,'assess_completed.py'),('Lung_case',P/'phase91_lung_patient_validation_20261015','review_lung_results.py')]
watch={};results=[]
def record():
 value={'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'waiting':[v[0] for v in watch.values()],'completed_events':results,'formal_patient_admission_automatic':False}
 f=BASE/'resources/DETACHED_REVIEW_EVENTS.json';tmp=f.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
for name,o,script in configs:
 pid=json.loads((o/'LAUNCH.json').read_text(encoding='utf-8'))['pid']
 try:fd=os.pidfd_open(pid);watch[fd]=(name,o,script)
 except ProcessLookupError:watch[-len(watch)-1]=(name,o,script)
record()
while watch:
 ready=[fd for fd in watch if fd<0] or select.select(list(watch),[],[])[0]
 for fd in ready:
  name,o,script=watch.pop(fd)
  if fd>=0:os.close(fd)
  result={'source':name,'observed_exit_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  try:
   done=json.loads((o/'JOB_FINISHED.json').read_text(encoding='utf-8'));result['job_exit_code']=done.get('exit_code')
   if done.get('exit_code')==0:
    log=o/'logs'/('detached_'+script+'.log')
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(o/'code'/script)],stdout=f,stderr=subprocess.STDOUT)
    result.update(review_exit_code=r.returncode,review_log=str(log),scientific_admission_pending=True)
   else:result['status']='analysis_failed_intervention_required'
  except Exception:result.update(status='event_review_error',traceback=traceback.format_exc())
  results.append(result);record()
record()
