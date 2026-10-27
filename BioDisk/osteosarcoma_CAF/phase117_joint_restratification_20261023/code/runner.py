from pathlib import Path
import os,sys,json,subprocess,datetime,hashlib,fcntl,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase117_joint_restratification_20261023')
for d in ['logs','receipts','resources']:(O/d).mkdir(parents=True,exist_ok=True)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python'
names={'01':'01_freeze.py','02':'02_structure.py','03':'03_correspondence.py','04':'04_function.py','05':'05_holdout.py','06':'06_acceptance_report.py'}
markers={'01':'FREEZE_COMPLETE.json','02':'STRUCTURE_COMPLETE.json','03':'CORRESPONDENCE_COMPLETE.json','04':'FUNCTION_COMPLETE.json','05':'HOLDOUT_COMPLETE.json','06':'FINAL_ACCEPTANCE.json'}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(v):
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(v,ensure_ascii=False)+'\n')
 print(json.dumps(v,ensure_ascii=False),flush=True)
stages=sys.argv[1:] or ['01','02','03','04','05','06'];final=0
for n in stages:
 script=O/'code'/names[n];assert script.exists(),script
 receipt=O/'receipts'/(n+'.json');dependencies=[script,O/'code/caf117_common.py']+([O/'code/caf117_function_utils.py'] if n in ['04','05'] else [])
 codehash={str(p):sha(p) for p in dependencies}
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old.get('exit_code')==0 and old.get('code_hashes')==codehash and (O/markers[n]).exists():emit(dict(utc=now(),stage=n,status='reused_verified_checkpoint'));continue
 attempt=len(list((O/'receipts').glob(n+'_attempt*.json')))+1;start=now();log=O/'logs'/(n+'_attempt'+str(attempt)+'.log');env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='4'
 emit(dict(utc=start,stage=n,status='native_started',script=str(script),attempt=attempt))
 with log.open('w',encoding='utf-8') as f:
  p=subprocess.Popen([PY,'-u',str(script)],cwd=O/'code',env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  (O/'active_child.json').write_text(json.dumps(dict(stage=n,pid=p.pid,started_at=start,log=str(log))),encoding='utf-8')
  rc=p.wait()
 good=rc==0 and (O/markers[n]).exists();rec=dict(stage=n,attempt=attempt,started_at=start,finished_at=now(),exit_code=rc,required_marker=markers[n],marker_exists=(O/markers[n]).exists(),code_hashes=codehash,log=str(log),status='completed' if good else 'failed')
 (O/'logs'/(n+'.exit_code')).write_text(str(rc)+'\n',encoding='utf-8');receipt.write_text(json.dumps(rec,indent=2),encoding='utf-8');(O/'receipts'/(n+'_attempt'+str(attempt)+'.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8');emit(rec)
 if not good:final=rc or 1;break
(O/'RUN_FINISHED.json').write_text(json.dumps(dict(utc=now(),requested_stages=stages,exit_code=final,all_requested_stages_finished=final==0)),encoding='utf-8');sys.exit(final)
