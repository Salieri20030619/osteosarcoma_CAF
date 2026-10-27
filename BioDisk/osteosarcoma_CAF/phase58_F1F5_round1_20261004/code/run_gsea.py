from pathlib import Path
import subprocess
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase58_F1F5_round1_20261004')
f=(p/'logs/gsea_pipeline.log').open('w')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/gsea_pipeline.py')],stdout=f,stderr=subprocess.STDOUT)
(p/'gsea_pipeline.exit_code').write_text(str(r.returncode))
