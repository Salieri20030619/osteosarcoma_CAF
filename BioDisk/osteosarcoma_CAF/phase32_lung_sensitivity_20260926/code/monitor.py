from pathlib import Path
import os,json,select,subprocess,datetime,fcntl,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase32_lung_sensitivity_20260926');PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
lock=(O/'monitor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
jobs={'hlca':(O/'hlca.pid',O/'HLCA_STATUS.json'),'sensitivity':(O/'sensitivity_fixed_features/sensitivity.pid',O/'sensitivity_fixed_features/SENSITIVITY_STATUS.json')}
def now():return datetime.datetime.now().astimezone().isoformat()
results={}
def write():
 stages={}
 for kind,(pidf,sf) in jobs.items():
  try:stages[kind]=json.loads(sf.read_text())
  except Exception:stages[kind]={'status':'UNREADABLE'}
 d={'at':now(),'status':'RUNNING' if len(results)<len(jobs) else 'ALL_QUEUED_TASKS_REVIEWED','completion_wait':'native pidfd events, 15-minute resource snapshots','stages':stages,'acceptance':results,'memory_available_kB':int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))}
 tmp=O/'MONITOR_STATUS.tmp';tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2));tmp.replace(O/'MONITOR_STATUS.json')
def validate(kind):
 with (O/'logs'/f'{kind}_acceptance.log').open('x') as log:
  p=subprocess.Popen([PY,str(O/'code/accept_outputs.py'),kind],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT);rc=p.wait()
 results[kind]={'validation_exit_code':rc,'accepted':rc==0,'at':now()}
fds={}
for kind,(pidf,sf) in jobs.items():
 try:fds[os.pidfd_open(int(pidf.read_text()))]=kind
 except ProcessLookupError:validate(kind)
write()
while fds:
 ready=select.select(list(fds),[],[],900)[0]
 for fd in ready:
  kind=fds.pop(fd);os.close(fd);validate(kind)
 write()
(O/'monitor.exit_code').write_text('0' if all(r['accepted'] for r in results.values()) else '1')
