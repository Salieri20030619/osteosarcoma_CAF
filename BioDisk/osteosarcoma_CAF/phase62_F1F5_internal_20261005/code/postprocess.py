from pathlib import Path
import subprocess,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005')
assert (O/'ANALYSIS_COMPLETE').exists()
for stage in ['summarize','review_function','build','verify']:
 with (O/'logs'/f'{stage}.log').open('w') as f:r=subprocess.run([sys.executable,'-u',str(O/'code'/f'{stage}.py')],stdout=f,stderr=subprocess.STDOUT)
 (O/f'{stage}.exit_code').write_text(str(r.returncode));print(stage,r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
