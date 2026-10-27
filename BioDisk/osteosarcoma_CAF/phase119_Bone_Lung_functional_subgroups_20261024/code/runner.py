from pathlib import Path
import subprocess,sys,json,datetime,hashlib,fcntl,os
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase119_Bone_Lung_functional_subgroups_20261024');lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
stages={'02':('02_structure.py','STRUCTURE_COMPLETE.json'),'03':('03_function.py','FUNCTION_COMPLETE.json'),'04':('04_holdout.py','HOLDOUT_COMPLETE.json'),'05':('05_verify_report.py','ACCEPTANCE.json')}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(x):
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(x,ensure_ascii=False)+'\n')
 print(json.dumps(x),flush=True)
rc=0
for key in sys.argv[1:]:
 name,marker=stages[key];p=O/'code'/name;files=[p,O/'code/caf119_common.py']
 if key in ['02','04']:files.append(O/'code/caf119_hierarchy.py')
 if key in ['03','04']:files.append(O/'code/caf119_statistics.py')
 hashes={str(f):sha(f) for f in files};receipt=O/'receipts'/(key+'.json')
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old.get('exit_code')==0 and old.get('code_hashes')==hashes and (O/marker).exists():record(dict(utc=now(),stage=key,status='reuse_completed'));continue
 attempt=len(list((O/'receipts').glob(key+'_attempt*.json')))+1;start=now();log=O/'logs'/(key+'_attempt'+str(attempt)+'.log');record(dict(utc=start,stage=key,status='started',attempt=attempt))
 with log.open('w',encoding='utf-8') as f:
  child=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python','-u',str(p)],cwd=O/'code',stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  (O/'ACTIVE_PROCESS.json').write_text(json.dumps(dict(stage=key,pid=child.pid,log=str(log),started_at=start)),encoding='utf-8');rc=child.wait()
 rec=dict(started_at=start,finished_at=now(),stage=key,attempt=attempt,exit_code=rc,marker_exists=(O/marker).exists(),code_hashes=hashes,log=str(log));receipt.write_text(json.dumps(rec,indent=2),encoding='utf-8');(O/'receipts'/(key+'_attempt'+str(attempt)+'.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8');(O/'logs'/(key+'.exit_code')).write_text(str(rc),encoding='utf-8');record(rec)
 if rc or not (O/marker).exists():rc=rc or 1;break
(O/'RUN_FINISHED.json').write_text(json.dumps(dict(finished_at=now(),requested=sys.argv[1:],exit_code=rc)),encoding='utf-8');sys.exit(rc)
