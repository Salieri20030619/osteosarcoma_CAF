from pathlib import Path
import subprocess,sys,json,datetime,pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase65_spatial_structure_20261005')
p=O/'code/verify.py';s=p.read_text().replace('low_memory=False)','low_memory=False,float_precision="round_trip")');p.write_text(s)
with (O/'logs/verify.log').open('w') as f:r=subprocess.run([sys.executable,'-u',str(p)],stdout=f,stderr=subprocess.STDOUT)
print('VERIFY_EXIT',r.returncode,(O/'logs/verify.log').read_text()[-4500:],flush=True)
(O/'QA/VERIFY_RUN_STATUS.json').write_text(json.dumps({'exit_status':r.returncode,'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reader':'round_trip to retain exact float ties'},indent=2))
mf=pd.read_csv(O/'resources/SPATIAL_INPUT_MANIFEST.tsv',sep='\t')
p=Path(mf.image_path.iloc[0]);print('IMAGE_PATH',p);print('SIBLINGS',[x.name for x in p.parent.iterdir()][:40],flush=True)
sys.exit(r.returncode)
