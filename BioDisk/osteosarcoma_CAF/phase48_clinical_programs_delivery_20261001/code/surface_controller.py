from pathlib import Path
import subprocess,sys,traceback
Q=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001')
if '--start' in sys.argv:
 f=(Q/'logs/surface_controller.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'surface_controller.pid').write_text(str(p.pid));print('SURFACE_STARTED',p.pid)
else:
 rc=1
 try:
  with (Q/'logs/surface.log').open('w') as f:p=subprocess.Popen([sys.executable,str(Q/'code/surface.py')],stdout=f,stderr=subprocess.STDOUT);(Q/'surface.pid').write_text(str(p.pid));rc=p.wait()
  (Q/'surface.exit_code').write_text(str(rc));assert rc==0 and (Q/'SURFACE_COMPLETE').read_text().strip()=='SUCCESS'
 except Exception:traceback.print_exc();rc=1
 (Q/'surface_controller.exit_code').write_text(str(rc));print('SURFACE_CONTROLLER_EXIT',rc)
