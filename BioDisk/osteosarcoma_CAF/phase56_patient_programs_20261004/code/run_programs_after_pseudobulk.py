from pathlib import Path
import os,select,subprocess,sys,json
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase56_patient_programs_20261004')
try:
 if not (p/'pseudobulk.exit_code').exists():
  pid=int((p/'pseudobulk.pid').read_text())
  try:
   fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 assert (p/'pseudobulk.exit_code').read_text().strip()=='0' and (p/'PSEUDOBULK_COMPLETE').is_file()
 f=(p/'logs/programs.log').open('w')
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/patient_program_statistics.py')],stdout=f,stderr=subprocess.STDOUT)
 (p/'programs.exit_code').write_text(str(r.returncode))
except Exception as e:
 (p/'programs_dependency_error.txt').write_text(repr(e));raise
