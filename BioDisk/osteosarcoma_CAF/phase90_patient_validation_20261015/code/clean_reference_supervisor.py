from pathlib import Path
import subprocess,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015')
with (O/'logs/NCC__BS7_clean_reference.log').open('w',encoding='utf-8') as f:
 r=subprocess.run(['bash',str(O/'code/run_native_clean.sh'),'NCC__BS7'],stdout=f,stderr=subprocess.STDOUT)
v={'exit_code':r.returncode,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completion_marker':(O/'cnv_clean/infercnv/NCC__BS7/COMPLETE').exists()}
(O/'resources/BS7_CLEAN_REFERENCE_EXIT.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
