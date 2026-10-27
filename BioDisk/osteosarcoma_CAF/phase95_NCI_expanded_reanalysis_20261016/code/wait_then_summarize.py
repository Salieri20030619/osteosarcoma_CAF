from pathlib import Path
import os,select,subprocess,json,fcntl,datetime
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase95_NCI_expanded_reanalysis_20261016')
f=(p/'completion_callback.lock').open('a');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 fd=os.pidfd_open(2324187);select.select([fd],[],[]);os.close(fd)
except ProcessLookupError:pass
if not (p/'COMPUTATION_FINISHED.json').exists():
 (p/'resources/CALLBACK_ERROR.json').write_text(json.dumps({'error':'Supervisor exited without computation completion receipt','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}),encoding='utf-8')
 raise SystemExit(1)
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(p/'code/nci_summarize.py')])
(p/'logs/summary.exit_code').write_text(str(r.returncode),encoding='utf-8')
