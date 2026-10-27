from pathlib import Path
import subprocess,os
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase59_four_functional_axes_20261004')
(O/'logs').mkdir(parents=True,exist_ok=True)
with (O/'logs/analysis.log').open('w') as f:
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/analyze.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'analysis.exit_code').write_text(str(r.returncode))
print('ANALYSIS_EXIT',r.returncode,flush=True)
raise SystemExit(r.returncode)
