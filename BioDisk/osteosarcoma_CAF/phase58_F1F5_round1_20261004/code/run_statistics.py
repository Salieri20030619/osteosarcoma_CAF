from pathlib import Path
import subprocess
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase58_F1F5_round1_20261004')
f=(p/'logs/statistics.log').open('w')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/statistics.py')],stdout=f,stderr=subprocess.STDOUT)
(p/'statistics.exit_code').write_text(str(r.returncode))
