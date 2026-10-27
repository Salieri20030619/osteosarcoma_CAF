from pathlib import Path
import subprocess
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005');(O/'logs').mkdir(parents=True,exist_ok=True)
for stage in ['prepare','gsea']:
 with (O/'logs'/f'{stage}.log').open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code'/f'{stage}.py')],stdout=f,stderr=subprocess.STDOUT)
 (O/f'{stage}.exit_code').write_text(str(r.returncode));print(stage,r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
(O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n')
