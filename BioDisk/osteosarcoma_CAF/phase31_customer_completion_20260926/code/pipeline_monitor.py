from pathlib import Path
import os,json,datetime,select,subprocess,traceback,fcntl
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase31_customer_completion_20260926')
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
def now():return datetime.datetime.now().astimezone().isoformat()
def jread(f):
 try:return json.loads(f.read_text())
 except Exception:return None
def snapshot(state):
 d={'at':now(),'status':state,'inspection_interval_seconds':900,'waiting':'native pidfd process completion events; periodic resource snapshots every 15 minutes','stages':{f.name:jread(f) for f in O.glob('*_STATUS.json') if f.name!='PIPELINE_STATUS.json'},'memory_available_kB':int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))}
 f=O/'PIPELINE_STATUS.json';t=f.with_suffix('.tmp');t.write_text(json.dumps(d,ensure_ascii=False,indent=2));t.replace(f)
 with (O/'logs/pipeline_events.jsonl').open('a') as log:log.write(json.dumps({'at':now(),'status':state})+chr(10))
def wait_pid(pid):
 try:fd=os.pidfd_open(pid)
 except ProcessLookupError:return
 try:
  while not select.select([fd],[],[],900)[0]:snapshot('WAITING_FOR_NATIVE_METHODS')
 finally:os.close(fd)
lock=(O/'pipeline.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 snapshot('WAITING_FOR_NATIVE_METHODS')
 wait_pid(int((O/'methods.pid').read_text()))
 # Process exit alone is insufficient. Preserve missing markers and upstream errors as evidence.
 verified=(O/'METHODS_FINISHED').exists() and (O/'methods.exit_code').exists() and (O/'delivery/METHODS_QA.json').exists()
 if not verified:
  snapshot('STOPPED_UPSTREAM_COMPLETION_UNVERIFIED');raise RuntimeError('Methods controller disappeared without complete exit status/markers/QA; no success inferred')
 f=O/'code/ambient_comparison.py'
 with (O/'logs/ambient_comparison.log').open('x') as log:
  p=subprocess.Popen([PY,'-u',str(f)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
  (O/'ambient.pid').write_text(str(p.pid));snapshot('AMBIENT_COMPARISON_RUNNING')
  try:fd=os.pidfd_open(p.pid)
  except ProcessLookupError:fd=None
  if fd is not None:
   while not select.select([fd],[],[],900)[0]:snapshot('AMBIENT_COMPARISON_RUNNING')
   os.close(fd)
  rc=p.wait()
 if rc or not (O/'ambient_comparison/COMPLETE').exists():raise RuntimeError('Ambient comparison failed or missing completion marker')
 snapshot('QUEUED_QUALITY_COMPARISON_FINISHED_REVIEW_REQUIRED');(O/'pipeline.exit_code').write_text('0\n')
except Exception:
 traceback.print_exc();snapshot('PIPELINE_FAILED_REVIEW_REQUIRED');(O/'pipeline.exit_code').write_text('1\n');raise
