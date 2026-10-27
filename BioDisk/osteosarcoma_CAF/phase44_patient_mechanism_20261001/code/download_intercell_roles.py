from pathlib import Path
import os,sys,json,requests,hashlib,datetime,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase44_patient_mechanism_20261001';R=O/'resources/intercell';R.mkdir(exist_ok=True)
if '--start' in sys.argv:
 assert not (O/'intercell_download.pid').exists(),'intercell download already launched'
 env=os.environ.copy()
 for k in list(env):
  if k.lower().endswith('_proxy'):env.pop(k,None)
 env['NO_PROXY']='*'
 with (O/'logs/intercell_download.log').open('x') as log:
  p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',__file__],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
 (O/'intercell_download.pid').write_text(str(p.pid)+'\n');print('DIRECT_INTERCELL_DOWNLOAD_STARTED',p.pid);raise SystemExit(0)
for k in list(os.environ):
 if k.lower().endswith('_proxy'):os.environ.pop(k,None)
s=requests.Session();s.trust_env=False;out=[];rc=0
for name,url in [('ligand_roles.tsv','https://omnipathdb.org/intercell?scope=generic&categories=ligand&causality=trans&license=commercial&format=tsv'),('receptor_roles.tsv','https://omnipathdb.org/intercell?scope=generic&categories=receptor&causality=rec&license=commercial&format=tsv')]:
 entry={'file':name,'url':url,'state':'downloading','proxy_used':False};out.append(entry)
 try:
  f=R/name;part=f.with_suffix('.part')
  with s.get(url,stream=True,timeout=(20,90)) as r:
   r.raise_for_status()
   with part.open('wb') as h:
    for b in r.iter_content(1048576):h.write(b)
  lines=part.read_text().splitlines();assert len(lines)>1 and len(lines[0].split('\t'))>2 and not lines[0].startswith('<'),'invalid role response';part.replace(f)
  entry.update(state='verified',bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),rows=len(lines)-1,columns=lines[0].split('\t'))
 except Exception as e:entry.update(state='failed',error=repr(e));rc=1
 (O/'INTERCELL_DOWNLOAD_STATUS.json').write_text(json.dumps({'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'resources':out,'state':'downloading'},indent=2));print(entry,flush=True)
(O/'INTERCELL_DOWNLOAD_STATUS.json').write_text(json.dumps({'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'resources':out,'state':'complete' if rc==0 else 'failed'},indent=2));(O/'intercell_download.exit_code').write_text(str(rc)+'\n')
if rc==0:(O/'INTERCELL_DOWNLOAD_COMPLETE').write_text('SUCCESS\n')
raise SystemExit(rc)