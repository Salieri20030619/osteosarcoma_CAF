from pathlib import Path
import sys,subprocess,json,traceback,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');Q=P/'phase48_clinical_programs_delivery_20261001';O=P/'phase46_spatial_functional_20261001'
if '--start' in sys.argv:
 f=(Q/'logs/controller.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'controller.pid').write_text(str(p.pid));print('CLINICAL_PROGRAMS_STARTED',p.pid)
else:
 rc=1
 try:
  for name,root,script in [('clinical_programs',Q,Q/'code/clinical_programs.py'),('spatial_sensitivity_repair',O,O/'code/sensitivity.py')]:
   (Q/'STATUS.json').write_text(json.dumps({'state':'running','step':name,'at':datetime.datetime.now().astimezone().isoformat()}))
   with (root/'logs'/f'{name}.log').open('w') as f:p=subprocess.Popen([sys.executable,str(script)],stdout=f,stderr=subprocess.STDOUT);(root/f'{name}.pid').write_text(str(p.pid));r=p.wait()
   (root/f'{name}.exit_code').write_text(str(r));assert r==0,(name,r)
  (Q/'STATUS.json').write_text('{"state":"computation_complete_QA_pending"}');rc=0
 except Exception as e:traceback.print_exc();(Q/'STATUS.json').write_text(json.dumps({'state':'failed','error':repr(e)}))
 (Q/'controller.exit_code').write_text(str(rc));print('CONTROLLER_EXIT',rc,flush=True)
