from pathlib import Path
import os,subprocess,json,datetime,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001')
if '--start' in sys.argv:
 f=(O/'logs/controller.log').open('a');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(O/'controller.pid').write_text(str(p.pid));print('STARTED',p.pid)
else:
 rc=1
 try:
  with (O/'logs/spatial.log').open('w') as f:p=subprocess.Popen([sys.executable,str(O/'code/spatial.py')],stdout=f,stderr=subprocess.STDOUT);(O/'spatial.pid').write_text(str(p.pid));rc=p.wait()
  (O/'spatial.exit_code').write_text(str(rc));assert rc==0 and (O/'SPATIAL_COMPLETE').read_text().strip()=='SUCCESS'
 except Exception as e:print(repr(e),flush=True)
 (O/'controller.exit_code').write_text(str(rc));print('CONTROLLER_FINISHED',rc,flush=True)
