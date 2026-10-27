from pathlib import Path
import json,datetime,shutil,subprocess,sys,concurrent.futures
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase65_spatial_structure_20261005')
p=O/'code/spatial_context.py';s=p.read_text()
s=s.replace("assert np.array_equal(lib,t.library_UMI) and (lib>0).all()","assert np.array_equal(lib,t.library_UMI)\n t['valid_library']=lib>0")
s=s.replace("/lib[:,None]*1e4","/np.where(lib>0,lib,np.nan)[:,None]*1e4")
p.write_text(s)
p=O/'code/discovery_robustness.py';s=p.read_text()
s=s.replace("'effect':float(dv[have].mean()) if ok else np.nan","'effect':(0.0 if abs(float(dv[have].mean()))<1e-12 else float(dv[have].mean())) if ok else np.nan")
p.write_text(s)
(O/'resources/NUMERICAL_QC_AMENDMENT.json').write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spatial':'MR4 contains one zero-library spot; barcode alignment and full libraries exactly agree with frozen source. Retain spot in exported roster; undefined normalized scores NaN; omit only from correlations, not patient count. Same rule applied to every sample.','discovery':'Absolute mean paired log2 effect below 1e-12 snapped to zero, avoiding phantom positive counts and invalid all-zero inference from numerical residue.','reason':'Input QC and numerical QC, no direction/P-based exclusion.'},indent=2))
shutil.copy2(O/'RUN_STATUS.json',O/'logs/INITIAL_RUN_STATUS_FAILED.json')
shutil.copy2(O/'logs/spatial_context.py.log',O/'logs/spatial_context_INITIAL_FAILED.log')
def run(name):
 log=O/'logs'/f'{name}.log'
 with log.open('w') as f:r=subprocess.run([sys.executable,'-u',str(O/'code'/name)],stdout=f,stderr=subprocess.STDOUT)
 print(name,r.returncode,log.read_text()[-2000:],flush=True)
 return {'script':name,'exit_status':r.returncode,'log':str(log)}
with concurrent.futures.ThreadPoolExecutor(2) as ex:steps=list(ex.map(run,['discovery_robustness.py','spatial_context.py']))
st={'success':all(x['exit_status']==0 for x in steps),'steps':steps,'finished_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'RUN_STATUS.json').write_text(json.dumps(st,indent=2));sys.exit(0 if st['success'] else 1)
