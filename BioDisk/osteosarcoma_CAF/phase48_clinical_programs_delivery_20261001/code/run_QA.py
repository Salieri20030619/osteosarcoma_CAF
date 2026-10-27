from pathlib import Path
import sys,subprocess,json,traceback
Q=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001');R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
if '--start' in sys.argv:
 f=(Q/'logs/QA_controller.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'QA_controller.pid').write_text(str(p.pid));print('QA_STARTED',p.pid)
else:
 rc=1
 try:
  assert (Q/'controller.exit_code').read_text().strip()=='0'
  with (Q/'logs/independent_R.log').open('w') as f:p=subprocess.Popen([R,'--vanilla',str(Q/'code/clinical_QA.R')],stdout=f,stderr=subprocess.STDOUT);(Q/'independent_R.pid').write_text(str(p.pid));rc=p.wait()
  (Q/'independent_R.exit_code').write_text(str(rc));assert rc==0
  assert (Q/'INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS';(Q/'PIPELINE_COMPLETE').write_text('SUCCESS');(Q/'STATUS.json').write_text('{"state":"complete","QA_accepted":true}')
 except Exception:traceback.print_exc();rc=1
 (Q/'QA_controller.exit_code').write_text(str(rc));print('QA_CONTROLLER_EXIT',rc,flush=True)
