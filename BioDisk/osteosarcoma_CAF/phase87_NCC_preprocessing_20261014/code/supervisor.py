from pathlib import Path
import subprocess, json, datetime
O=Path(__file__).resolve().parents[1]
with (O/'logs/pipeline.log').open('w',encoding='utf-8') as f:
    p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python','-u',str(O/'code/process_ncc.py')],stdout=f,stderr=subprocess.STDOUT)
    (O/'worker.pid').write_text(str(p.pid)+'\n',encoding='utf-8')
    rc=p.wait()
(O/'logs/pipeline.exit_code').write_text(str(rc)+'\n',encoding='utf-8')
ok=rc==0 and (O/'PREPROCESSING_COMPLETE.json').exists()
(O/'JOB_FINISHED.json').write_text(json.dumps({'exit_code':rc,'completion_marker_present':(O/'PREPROCESSING_COMPLETE.json').exists(),'success':ok,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n',encoding='utf-8')
