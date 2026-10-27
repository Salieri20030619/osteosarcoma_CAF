from pathlib import Path
import subprocess,os,json,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase41_full_transcriptome_setup_20261001'
for sub in ['resources','logs','code','tables','objects','delivery']:(O/sub).mkdir(parents=True,exist_ok=True)
f=O/'code/download_resources.py'
assert not (O/'download.pid').exists(), 'Download already launched; inspect status instead'
env=os.environ.copy()
for key in list(env):
 if key.lower().endswith('_proxy'):env.pop(key,None)
env['NO_PROXY']='*'
with (O/'logs/download.log').open('x') as log:
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(f)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
(O/'download.pid').write_text(str(p.pid)+'\n')
(O/'DOWNLOAD_LAUNCH.json').write_text(json.dumps({'pid':p.pid,'proxy_used':False,'ssh_alias':'sgs10','at':datetime.datetime.now().astimezone().isoformat(),'script':str(f)},indent=2))
print('DIRECT_BACKGROUND_DOWNLOAD_STARTED',p.pid)