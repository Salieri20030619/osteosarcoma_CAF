from pathlib import Path
import subprocess,sys,json,datetime,concurrent.futures,shutil
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase65_spatial_structure_20261005')
shutil.copy2(O.parent/'tmp_phase65_freeze.py',O/'code/freeze.py')
def run(name):
 p=O/'logs'/f'{name}.log'
 with p.open('w') as f:r=subprocess.run([sys.executable,'-u',str(O/'code'/name)],stdout=f,stderr=subprocess.STDOUT)
 result={'script':name,'exit_status':r.returncode,'log':str(p)}
 print('FINISHED',result,flush=True);print(p.read_text()[-2000:],flush=True)
 return result
status={'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
with concurrent.futures.ThreadPoolExecutor(2) as ex:status['steps']=list(ex.map(run,['discovery_robustness.py','spatial_context.py']))
status['success']=all(x['exit_status']==0 for x in status['steps']);status['finished_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(O/'RUN_STATUS.json').write_text(json.dumps(status,indent=2));sys.exit(0 if status['success'] else 1)
