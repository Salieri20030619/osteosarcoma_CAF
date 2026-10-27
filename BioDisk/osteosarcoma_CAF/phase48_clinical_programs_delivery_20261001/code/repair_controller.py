from pathlib import Path
import sys,subprocess,json,traceback
Q=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001');R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
if '--start' in sys.argv:
 f=(Q/'logs/repair_controller.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'repair_controller.pid').write_text(str(p.pid));print('REPAIR_STARTED',p.pid)
else:
 rc=1
 try:
  for name,cmd in [('clinical_programs_repaired',[sys.executable,str(Q/'code/clinical_programs.py')]),('independent_R_repaired',[R,'--vanilla',str(Q/'code/clinical_QA.R')])]:
   with (Q/'logs'/f'{name}.log').open('w') as f:p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT);(Q/f'{name}.pid').write_text(str(p.pid));r=p.wait()
   (Q/f'{name}.exit_code').write_text(str(r));assert r==0,(name,r)
  assert (Q/'INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS';(Q/'PIPELINE_COMPLETE').write_text('SUCCESS');(Q/'STATUS.json').write_text('{"state":"complete","HC3_unit_leverage_excluded":true,"QA_accepted":true}');rc=0
 except Exception:traceback.print_exc()
 (Q/'repair_controller.exit_code').write_text(str(rc));print('REPAIR_EXIT',rc,flush=True)
