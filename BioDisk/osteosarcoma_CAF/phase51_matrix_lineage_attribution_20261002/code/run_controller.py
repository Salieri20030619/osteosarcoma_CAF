from pathlib import Path
import subprocess,json
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase51_matrix_lineage_attribution_20261002')
cmd=['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(r/'code/analysis.py')]
(r/'resources/RUN_COMMAND.json').write_text(json.dumps(cmd))
with (r/'logs/analysis.log').open('w') as log:
 proc=subprocess.Popen(cmd,cwd=r,stdout=log,stderr=subprocess.STDOUT);(r/'analysis.pid').write_text(str(proc.pid));rc=proc.wait()
(r/'analysis.exit_code').write_text(str(rc))
print(json.dumps({'exit_code':rc,'marker':(r/'ANALYSIS_COMPLETE').exists(),'log_tail':(r/'logs/analysis.log').read_text()[-5000:]}))