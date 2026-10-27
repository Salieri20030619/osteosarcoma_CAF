from pathlib import Path
import os,ctypes,subprocess,json,traceback,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase49_audit_revision_20261002');target=R/'pathway_revision.exit_code'
rc=1
try:
 if not target.exists():
  libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(0);assert fd>=0;assert libc.inotify_add_watch(fd,os.fsencode(R),8|128)>=0
  while not target.exists():os.read(fd,65536)
  os.close(fd)
 assert target.read_text().strip()=='0' and (R/'PATHWAY_REVISION_COMPLETE').read_text().strip()=='SUCCESS'
 with (R/'logs/pathway_independent_QA.log').open('w') as f:p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(R/'code/verify_pathways.py')],stdout=f,stderr=subprocess.STDOUT)
 (R/'pathway_independent_QA.exit_code').write_text(str(p.returncode));assert p.returncode==0
 (R/'COMPUTATION_VALIDATED').write_text('SUCCESS');rc=0
except Exception:traceback.print_exc()
(R/'compute_validation_controller.exit_code').write_text(str(rc));print('COMPUTE_VALIDATION_EXIT',rc,flush=True)
