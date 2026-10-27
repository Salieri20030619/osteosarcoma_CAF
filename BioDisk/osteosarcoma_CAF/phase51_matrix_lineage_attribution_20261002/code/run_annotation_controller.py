from pathlib import Path
import subprocess,json
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase51_matrix_lineage_attribution_20261002')
with (r/'logs/annotation_CNV.log').open('w') as log:
 proc=subprocess.Popen(['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(r/'code/annotation_CNV_detail.py')],stdout=log,stderr=subprocess.STDOUT,cwd=r);(r/'annotation_CNV.pid').write_text(str(proc.pid));rc=proc.wait()
(r/'annotation_CNV.exit_code').write_text(str(rc))
print(json.dumps({'exit_code':rc,'marker':(r/'ANNOTATION_CNV_COMPLETE').exists(),'log_tail':(r/'logs/annotation_CNV.log').read_text()[-2000:]}))