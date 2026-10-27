from pathlib import Path
import subprocess,psutil
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase10_batch_gates_20260919')
if not (O/'EARLY_AMBIENT_COMPLETE').exists():
 jobs=[]
 for p in psutil.process_iter(['pid','cmdline']):
  cmd=p.info['cmdline'] or []
  if any(str(x).endswith('/complete_early_ambient.py') for x in cmd):jobs.append(p.info['pid'])
 assert jobs,'Missing ambient completion and no active producer'
 for pid in jobs:subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
assert (O/'EARLY_AMBIENT_COMPLETE').exists()
with open(O/'logs/unify_gates.log','w') as f:rc=subprocess.call(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/unify_gates.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'logs/unify_gates.exit_code').write_text(str(rc)+'\n');print('UNIFY_GATES_EXIT',rc,flush=True)
if rc:print((O/'logs/unify_gates.log').read_text()[-6000:])
raise SystemExit(rc)