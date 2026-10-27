from pathlib import Path
import os,select,subprocess,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase91_lung_patient_validation_20261015')
if not (O/'PREPARATION_COMPLETE.json').exists():
 parent=json.loads((O/'LAUNCH.json').read_text(encoding='utf-8'))['pid']
 children=subprocess.check_output(['ps','--ppid',str(parent),'-o','pid=,args='],text=True)
 found=False
 for line in children.splitlines():
  if 'prepare_lung.py' in line:
   pid=int(line.split()[0]);found=True
   try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
   except ProcessLookupError:pass
 assert found or (O/'PREPARATION_COMPLETE.json').exists(),'No active preparation process or completion marker'
assert (O/'PREPARATION_COMPLETE.json').exists(),'Preparation failed'
with (O/'logs/fingerprint.log').open('w',encoding='utf-8') as f:
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(O/'code/fingerprint.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'resources/FINGERPRINT_EXIT.json').write_text(json.dumps({'exit_code':r.returncode,'completion_marker':(O/'resources/FINGERPRINT_COMPLETE.json').exists(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n',encoding='utf-8')
