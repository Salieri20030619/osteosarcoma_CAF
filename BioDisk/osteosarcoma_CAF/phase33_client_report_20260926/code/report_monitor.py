from pathlib import Path
import os,select,json,subprocess,datetime,fcntl
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase33_client_report_20260926')
lock=(O/'report_monitor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now().astimezone().isoformat()
def snapshot(s):
 f=O/'REPORT_MONITOR_STATUS.json';t=f.with_suffix('.tmp');d={'at':now(),'status':s,'native_completion_wait':True,'snapshot_interval_seconds':900}
 if (O/'REFERENCE_STATUS.json').exists():d['reference']=json.loads((O/'REFERENCE_STATUS.json').read_text())
 t.write_text(json.dumps(d,indent=2));t.replace(f)
pid=int((O/'reference.pid').read_text());snapshot('WAITING_FOR_REFERENCE_COMPARISON')
try:
 fd=os.pidfd_open(pid)
 while not select.select([fd],[],[],900)[0]:snapshot('WAITING_FOR_REFERENCE_COMPARISON')
 os.close(fd)
except ProcessLookupError:pass
assert (O/'REFERENCE_FINISHED').exists() and (O/'reference.exit_code').exists()
with (O/'logs/reference_comparison.log').open('x') as log:
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/compare_reference.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
 rc=p.wait()
(O/'comparison.exit_code').write_text(str(rc));snapshot('REFERENCE_REVIEW_READY' if rc==0 else 'COMPARISON_FAILED')
