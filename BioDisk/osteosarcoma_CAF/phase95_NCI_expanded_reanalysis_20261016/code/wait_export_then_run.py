from pathlib import Path
import os,select,subprocess,json,datetime,fcntl
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase95_NCI_expanded_reanalysis_20261016')
lock=(p/'bootstrap.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 fd=os.pidfd_open(2189124);select.select([fd],[],[]);os.close(fd)
except ProcessLookupError:pass
r=json.loads((p/'resources/EXPORT_JOB_FINISHED.json').read_text(encoding='utf-8'));assert r['exit_code']==0 and r['completion_marker']
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(p/'code/runner.py')])
(p/'logs/runner.exit_code').write_text(str(r.returncode),encoding='utf-8')
