
from pathlib import Path
import subprocess,os,json,datetime
from concurrent.futures import ThreadPoolExecutor
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase33_client_report_20260926')
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
def run(var):
 d=O/'cnv_reference/SCPCS000610'/var;d.mkdir(parents=True,exist_ok=True)
 env=os.environ.copy()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
 with (O/'logs'/('early_SCPCS000610_'+var+'.log')).open('x') as log:
  p=subprocess.Popen([R,'--vanilla',str(O/'code/reference_sensitivity.R'),'SCPCS000610',var],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env)
  (d/'early.pid').write_text(str(p.pid));rc=p.wait()
 (d/'early.exit_code').write_text(str(rc));return {'variant':var,'exit_code':rc,'complete':(d/'COMPLETE').exists()}
with ThreadPoolExecutor(max_workers=2) as pool:
 outcomes=list(pool.map(run,['source_reference','lymphoid_reference']))
(O/'early_reference_outcomes.json').write_text(json.dumps(outcomes,indent=2))
