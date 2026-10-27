from pathlib import Path
import os,select,subprocess,json,datetime,fcntl
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase92_copykat_validation_20261015')
lock=(p/'patient_verification.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 fd=os.pidfd_open(1227619);select.select([fd],[],[]);os.close(fd)
except ProcessLookupError:pass
assert (p/'COMPUTATION_RECONCILED.json').exists(),'Monitor exited without reconciled completion'
with (p/'logs/patient_count_verification.log').open('w',encoding='utf-8') as log:
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(p/'code/verify_patient_counts.py')],stdout=log,stderr=subprocess.STDOUT)
record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':r.returncode,'completion_marker':(p/'patient_count_verification/VERIFIED_COUNTS.json').exists(),'scientific_admission_pending':True}
(p/'resources/PATIENT_COUNT_VERIFICATION_EXIT.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print(json.dumps(record),flush=True)
