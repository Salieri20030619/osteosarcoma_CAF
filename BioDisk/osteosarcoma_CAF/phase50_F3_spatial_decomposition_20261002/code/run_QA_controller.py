from pathlib import Path
import subprocess,json
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase50_F3_spatial_decomposition_20261002')
assert (p/'analysis.exit_code').read_text().strip()=='0'
(p/'code/parse_QA.R').write_text('parse(file="'+str(p/'code/independent_QA.R')+'")')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(p/'code/parse_QA.R')],capture_output=True,text=True)
(p/'logs/R_parse.log').write_text(r.stdout+r.stderr)
assert r.returncode==0,r.stderr
with (p/'logs/independent_R.log').open('w') as log:
 proc=subprocess.Popen(['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(p/'code/independent_QA.R')],stdout=log,stderr=subprocess.STDOUT,cwd=p)
 (p/'independent_R.pid').write_text(str(proc.pid));rc=proc.wait()
(p/'independent_R.exit_code').write_text(str(rc))
print(json.dumps({'exit_code':rc,'marker':(p/'INDEPENDENT_R_QA_COMPLETE').exists(),'log_tail':(p/'logs/independent_R.log').read_text()[-5000:]}))