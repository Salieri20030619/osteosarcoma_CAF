from pathlib import Path
import subprocess
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase58_F1F5_round1_20261004')
with (p/'logs/QA_controller_retry.log').open('w') as f: r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(p/'code/prepare_QA.py')],stdout=f,stderr=subprocess.STDOUT)
(p/'QA.exit_code').write_text(str(r.returncode))
