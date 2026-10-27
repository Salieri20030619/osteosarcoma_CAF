from pathlib import Path
import subprocess
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase56_patient_programs_20261004')
f=(p/'logs/pseudobulk.log').open('w')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/prepare_patient_programs.py')],stdout=f,stderr=subprocess.STDOUT)
(p/'pseudobulk.exit_code').write_text(str(r.returncode))
