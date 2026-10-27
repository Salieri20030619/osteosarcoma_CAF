from pathlib import Path
import subprocess,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015')
with (O/'logs/fingerprint.log').open('w',encoding='utf-8') as f:
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(O/'code/fingerprint.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'resources/FINGERPRINT_EXIT.json').write_text(json.dumps({'exit_code':r.returncode,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completion_marker':(O/'resources/FINGERPRINT_COMPLETE.json').exists()},indent=2)+'\n',encoding='utf-8')
