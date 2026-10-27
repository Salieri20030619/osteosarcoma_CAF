from pathlib import Path
import subprocess,sys,py_compile
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase63_matrix_resistance_20261005')
for stage in sys.argv[1:]:
 script=O/'code'/f'{stage}.py';py_compile.compile(str(script),doraise=True)
 with (O/'logs'/f'{stage}.log').open('w') as f:r=subprocess.run([sys.executable,'-u',str(script)],stdout=f,stderr=subprocess.STDOUT)
 (O/f'{stage}.exit_code').write_text(str(r.returncode));print(stage,r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
