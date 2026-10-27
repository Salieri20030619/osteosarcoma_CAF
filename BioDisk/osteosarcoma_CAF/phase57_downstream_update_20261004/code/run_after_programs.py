from pathlib import Path
import os,select,subprocess
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');b=r/'phase56_patient_programs_20261004';o=r/'phase57_downstream_update_20261004'
try:
 if not (b/'programs.exit_code').exists():
  try:
   fd=os.pidfd_open(int((b/'programs_controller.pid').read_text()));select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 assert (b/'programs.exit_code').read_text().strip()=='0' and (b/'PROGRAMS_COMPLETE').is_file()
 with (o/'logs/downstream.log').open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(o/'code/downstream_update.py')],stdout=f,stderr=subprocess.STDOUT)
 (o/'downstream.exit_code').write_text(str(r.returncode))
except Exception as e:
 (o/'dependency_error.txt').write_text(repr(e));raise
