from pathlib import Path
import os,select,json,subprocess,datetime,time,traceback,fcntl
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015')
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
def dump(v):
 v['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'SECONDARY_PROGRESS.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
lock=(O/'secondary.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
pid=json.loads((O/'resources/BS7_CLEAN_REFERENCE_LAUNCH.json').read_text(encoding='utf-8'))['pid']
dump({'status':'waiting_native_clean_reference_completion','wait_pid':pid,'pid':os.getpid()})
try:
 try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
 receipt=json.loads((O/'resources/BS7_CLEAN_REFERENCE_EXIT.json').read_text(encoding='utf-8'));assert receipt['exit_code']==0 and receipt['completion_marker'],'Clean-reference run needs intervention first'
 jobs=[]
 for sm in ['NCC__BS9','NCC__BS5','NCC__BS4','NCC__BS3']:
  d=O/'inputs'/sm
  for name,cmd,marker in [('doublets',[R,'--vanilla',str(O/'code/doublets_full.R'),str(d)],d/'DOUBLETS_COMPLETE'),('infercnv',['bash',str(O/'code/run_native.sh'),sm],O/'cnv/infercnv'/sm/'COMPLETE')]:
   if marker.exists():continue
   while int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))<96*1024*1024:
    dump({'status':'waiting_for_memory','sample':sm,'method':name,'jobs':jobs});time.sleep(900)
   dump({'status':'running','sample':sm,'method':name,'pid':os.getpid(),'jobs':jobs})
   with (O/'logs'/(sm+'_'+name+'.log')).open('w',encoding='utf-8') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   (O/'logs'/(sm+'_'+name+'.exit_code')).write_text(str(r.returncode)+'\n',encoding='utf-8')
   z={'sample':sm,'method':name,'exit_code':r.returncode,'completion_marker':marker.exists()};jobs.append(z)
   assert r.returncode==0 and marker.exists(),z
 dump({'status':'complete','exit_code':0,'jobs':jobs})
except Exception:
 dump({'status':'failed','exit_code':1,'traceback':traceback.format_exc()});raise
