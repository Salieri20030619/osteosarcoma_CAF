from pathlib import Path
import os,subprocess,json,datetime,concurrent.futures,traceback,fcntl,time
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase91_lung_patient_validation_20261015')
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript';sm='GSE234187__GSM7453693'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
def dump(n,v):
 v['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/n).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def run(name,cmd,marker):
 with (O/'logs'/(name+'.log')).open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 (O/'logs'/(name+'.exit_code')).write_text(str(r.returncode)+'\n',encoding='utf-8')
 return {'job':name,'exit_code':r.returncode,'completion_marker':marker.exists()}
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 dump('STATUS.json',{'status':'preparing_full_input_and_doublets','pid':os.getpid()})
 r=run('prepare',[PY,str(O/'code/prepare_lung.py')],O/'PREPARATION_COMPLETE.json');assert r['exit_code']==0 and r['completion_marker'],r
 while int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))<96*1024*1024:time.sleep(900)
 dump('STATUS.json',{'status':'native_CNV_and_ambient_running','pid':os.getpid()})
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  fs=[ex.submit(run,'infercnv',['bash',str(O/'code/run_native.sh'),sm],O/'cnv/infercnv'/sm/'COMPLETE'),ex.submit(run,'decontx',[R,'--vanilla',str(O/'code/lung_decontx.R'),str(O/'inputs'/sm)],O/'ambient/COMPLETE')]
  results=[f.result() for f in fs]
 ok=all(r['exit_code']==0 and r['completion_marker'] for r in results)
 dump('JOB_FINISHED.json',{'exit_code':0 if ok else 1,'jobs':results,'scientific_admission_pending':True,'new_formal_patients':0})
 dump('STATUS.json',{'status':'computation_complete_review_pending' if ok else 'tool_failure','new_formal_patients':0})
except Exception:
 dump('JOB_FINISHED.json',{'exit_code':1,'traceback':traceback.format_exc()});dump('STATUS.json',{'status':'failed'});raise
