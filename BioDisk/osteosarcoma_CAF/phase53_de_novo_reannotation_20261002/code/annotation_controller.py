from pathlib import Path
import subprocess,json,datetime,os,ctypes,concurrent.futures,threading
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002')
libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(os.O_CLOEXEC)
if fd<0:raise OSError(ctypes.get_errno(),"inotify_init1")
wd=libc.inotify_add_watch(fd,str(r).encode(),0x8|0x80|0x100)
if wd<0:raise OSError(ctypes.get_errno(),"inotify_add_watch")
while not (r/'structure_pipeline.exit_code').exists():os.read(fd,65536)
os.close(fd)
if (r/'structure_pipeline.exit_code').read_text().strip()!='0':
 (r/'annotation_pipeline.exit_code').write_text('1');raise SystemExit("Structure failure; annotation not started")
datasets=['SCPCP000017','SCPCP000018','SCPCP000023','GSE152048','GSE318840','GSE162454','GSE270231','GSE299023','GSE250015']
state={};lock=threading.Lock()
def record(ds,v):
 with lock:
  state[ds]=v;p=r/'ANNOTATION_PROGRESS.json';t=p.with_suffix('.tmp');t.write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'datasets':state},indent=2));t.replace(p)
def job(ds):
 marker=r/'objects'/ds/'REFERENCE_ANNOTATION_COMPLETE';exitfile=r/'logs'/f'annotation_{ds}.exit_code'
 if marker.exists() and exitfile.exists() and exitfile.read_text().strip()=='0':record(ds,{'state':'verified_completed_reuse'});return 0
 record(ds,{'state':'running'})
 with (r/'logs'/f'annotation_{ds}.log').open('w') as log:
  p=subprocess.Popen(['/usr/bin/time','-v','/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(r/'code/annotate_references.R'),ds],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=r)
  (r/'logs'/f'annotation_{ds}.pid').write_text(str(p.pid));rc=p.wait()
 exitfile.write_text(str(rc));good=rc==0 and marker.exists();record(ds,{'state':'complete' if good else 'failed','exit':rc});return int(not good)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(job,datasets))
rc=int(any(results));(r/'annotation_pipeline.exit_code').write_text(str(rc))
if rc==0:(r/'ANNOTATION_PIPELINE_COMPLETE').write_text('SUCCESS')
print(json.dumps({'exit':rc,'states':state}),flush=True)
