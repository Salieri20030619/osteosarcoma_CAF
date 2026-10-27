from pathlib import Path
import os,subprocess,json,datetime,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase43_full_transcriptome_discovery_20261001';PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
for sub in ['logs','delivery','resources']:(O/sub).mkdir(exist_ok=True)
assert not (O/'controller.pid').exists(),'Existing controller launch; inspect before retry'
env=os.environ.copy()
for key in list(env):
 if key.lower().endswith('_proxy'):env.pop(key,None)
env['NO_PROXY']='*'
(O/'resources/CODE_SHA256.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (O/'code').iterdir() if f.is_file()},indent=2))
with (O/'logs/controller.log').open('x') as log:
 p=subprocess.Popen([PY,'-u',str(O/'code/stageC_controller.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
(O/'controller.pid').write_text(str(p.pid)+'\n');(O/'delivery/LAUNCH.json').write_text(json.dumps({'pid':p.pid,'at':datetime.datetime.now().astimezone().isoformat(),'ssh_alias':'sgs10','proxy_used':False,'sequence':'C1 full genes -> independent R QA -> C2 ORA/GSEA','native_wait':True},indent=2));print('SEQUENTIAL_STAGE_C_STARTED',p.pid)