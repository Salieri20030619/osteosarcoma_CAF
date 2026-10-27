from pathlib import Path
import subprocess,sys,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
assert (O/'sources/GPL14951_platform_only.soft').exists(),'Complete GPL14951 annotation required'
status={'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'steps':[]}
for name in ['analyze_response.py','verify_response.py']:
 log=O/'logs'/f'{name}.log'
 with log.open('w') as f:r=subprocess.run([sys.executable,str(O/'code'/name)],stdout=f,stderr=subprocess.STDOUT)
 status['steps'].append({'script':name,'exit_status':r.returncode,'log':str(log)})
 print(name,'EXIT',r.returncode,flush=True);print(log.read_text()[-15000:],flush=True)
 if r.returncode:break
status['finished_at']=datetime.datetime.now(datetime.timezone.utc).isoformat();status['success']=len(status['steps'])==2 and all(s['exit_status']==0 for s in status['steps'])
(O/'FINAL_RUN_STATUS.json').write_text(json.dumps(status,indent=2));sys.exit(0 if status['success'] else 1)
