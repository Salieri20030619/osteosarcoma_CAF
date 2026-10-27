from pathlib import Path
import subprocess,json,datetime
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase52_blind_lineage_validation_20261002')
cmd=['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(r/'code/analysis.py')]
(r/'resources/RUN_COMMAND.json').write_text(json.dumps(cmd))
with (r/'logs/analysis.log').open('w') as log:
 p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=r)
 (r/'analysis.pid').write_text(str(p.pid))
 rc=p.wait()
(r/'analysis.exit_code').write_text(str(rc))
print(json.dumps({'exit_code':rc,'marker':(r/'ANALYSIS_COMPLETE').exists(),'log_tail':(r/'logs/analysis.log').read_text()[-4500:]}),flush=True)
