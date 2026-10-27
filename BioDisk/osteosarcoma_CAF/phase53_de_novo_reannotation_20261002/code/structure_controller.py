from pathlib import Path
import subprocess,json,datetime,concurrent.futures,threading,time,shutil,os
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002')
datasets=['SCPCP000017','SCPCP000018','SCPCP000023','GSE152048','GSE318840','GSE162454','GSE270231','GSE299023','GSE250015']
lock=threading.Lock();state={}
def record(ds,v):
 with lock:
  state[ds]=v;p=r/'STRUCTURE_PROGRESS.json';t=p.with_suffix('.tmp');t.write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'datasets':state},indent=2));t.replace(p)
def job(ds):
 marker=r/'objects'/ds/'STRUCTURE_COMPLETE';exitfile=r/'logs'/f'structure_{ds}.exit_code'
 if marker.exists() and exitfile.exists() and exitfile.read_text().strip()=='0':
  record(ds,{'state':'verified_completed_reuse'});return 0
 while True:
  available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
  if available>98*1024**3 and shutil.disk_usage(r).free>20*1024**3:break
  record(ds,{'state':'waiting_for_resources','available_GiB':available/1024**3});time.sleep(900)
 record(ds,{'state':'running'})
 with (r/'logs'/f'structure_{ds}.log').open('w') as log:
  p=subprocess.Popen(['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(r/'code/build_structure.py'),ds],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=r)
  (r/'logs'/f'structure_{ds}.pid').write_text(str(p.pid));rc=p.wait()
 exitfile.write_text(str(rc));good=rc==0 and marker.exists()
 record(ds,{'state':'complete' if good else 'failed','exit_code':rc,'marker':marker.exists()})
 return 0 if good else 1
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 results=list(pool.map(job,datasets))
rc=int(any(results));(r/'structure_pipeline.exit_code').write_text(str(rc))
if rc==0:(r/'STRUCTURE_PIPELINE_COMPLETE').write_text('SUCCESS')
print(json.dumps({'exit':rc,'states':state}),flush=True)
