from pathlib import Path
import subprocess,os,json
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase50_F3_spatial_decomposition_20261002')
with (p/'logs/analysis.log').open('w') as log:
 proc=subprocess.Popen(['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(p/'code/analysis.py')],stdout=log,stderr=subprocess.STDOUT,cwd=p)
 (p/'analysis.pid').write_text(str(proc.pid))
 rc=proc.wait()
(p/'analysis.exit_code').write_text(str(rc))
print(json.dumps({'exit_code':rc,'marker':(p/'ANALYSIS_COMPLETE').exists(),'log_tail':(p/'logs/analysis.log').read_text()[-6000:]}))