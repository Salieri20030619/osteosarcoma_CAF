from pathlib import Path
import subprocess,json,datetime,shutil,sys,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase120_patient_supported_CAF_review_20261026'
for s in ['code','logs','resources']: (O/s).mkdir(parents=True,exist_ok=True)
p=O/'code/01_audit.py';shutil.copyfile(R/'caf120_audit.py',p)
start=datetime.datetime.now(datetime.timezone.utc).isoformat()
with (O/'logs/01_audit.log').open('w',encoding='utf-8') as log:
 c=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python','-u',str(p)],stdout=log,stderr=subprocess.STDOUT)
 (O/'ACTIVE_PROCESS.json').write_text(json.dumps(dict(pid=c.pid,started_at=start,code_sha256=hashlib.sha256(p.read_bytes()).hexdigest())),encoding='utf-8')
 rc=c.wait()
rec=dict(exit_code=rc,started_at=start,finished_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),marker=(O/'AUDIT_COMPLETE.json').exists())
(O/'RUN_RECEIPT.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec),flush=True);sys.exit(rc)
