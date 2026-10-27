from pathlib import Path
import subprocess
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase54_raw_review_20261004')
f=(p/'logs/raw_review.log').open('w')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/raw_review.py')],stdout=f,stderr=subprocess.STDOUT)
(p/'raw_review.exit_code').write_text(str(r.returncode))
