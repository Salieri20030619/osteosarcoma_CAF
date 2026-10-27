from pathlib import Path
import subprocess,psutil,json,shutil,datetime,py_compile
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919'
py_compile.compile(str(O/'code/batch_models.py'),doraise=True)
if not (O/'UNIFIED_GATES_COMPLETE').exists():
 jobs=[]
 for p in psutil.process_iter(['pid','cmdline']):
  cmd=p.info['cmdline'] or []
  if any(str(x).endswith('/unify_gates.py') for x in cmd):jobs.append(p.info['pid'])
 assert jobs,'No completed unified gates or active producer'
 for pid in jobs:subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
assert (O/'UNIFIED_GATES_COMPLETE').exists()
statepath=P/'PROJECT_STATE.json';backup=P/'PROJECT_STATE.before_batch_gates_20260919.json'
if not backup.exists():shutil.copy2(statepath,backup)
s=json.loads(statepath.read_text());s.update(current_stage=O.name,status='BATCH_AND_GATE_REMEDIATION_RUNNING',current_status='BATCH_AND_GATE_REMEDIATION_RUNNING',dynamic_status_file=str(O/'STATUS.json'),next_action='Review unified gates and batch sensitivity models, preserve SR1004 exact target name');statepath.write_text(json.dumps(s,ensure_ascii=False,indent=2))
with open(O/'logs/batch_models.log','w') as f:rc=subprocess.call(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/batch_models.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'logs/batch_models.exit_code').write_text(str(rc)+'\n');print('BATCH_MODELS_EXIT',rc,flush=True)
if rc:print((O/'logs/batch_models.log').read_text()[-6000:])
else:print((O/'STATUS.json').read_text())
raise SystemExit(rc)