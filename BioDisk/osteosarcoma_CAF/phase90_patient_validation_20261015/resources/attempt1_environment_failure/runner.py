from pathlib import Path
import os,subprocess,json,datetime,traceback,time,fcntl
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase90_patient_validation_20261015'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
def dump(name,v):
 v['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 f=O/name;t=f.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');t.replace(f)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
jobs=[]
def run(name,args,marker):
 if marker.exists():return
 while int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))<96*1024*1024:
  dump('STATUS.json',{'status':'waiting_for_memory','next_job':name,'pid':os.getpid()});time.sleep(900)
 dump('STATUS.json',{'status':'running','job':name,'pid':os.getpid(),'new_formal_patients':0})
 with (O/'logs'/(name+'.log')).open('w',encoding='utf-8') as f:
  try:rc=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=16*3600).returncode
  except subprocess.TimeoutExpired:rc=124
 jobs.append({'job':name,'exit_code':rc,'completion_marker':marker.exists()});dump('JOB_PROGRESS.json',{'jobs':jobs})
 (O/'logs'/(name+'.exit_code')).write_text(str(rc)+'\n',encoding='utf-8')
 if rc or not marker.exists():raise RuntimeError(name+' failed; see log')
try:
 run('prepare',[PY,str(O/'code/prepare_full.py')],O/'PREPARATION_COMPLETE.json')
 import csv
 rows=list(csv.DictReader((O/'resources/manifest.tsv').open(encoding='utf-8'),delimiter='\t'))
 for r in rows:
  sm=r['sample'];d=O/'inputs'/sm
  run(sm+'_doublets',[R,'--vanilla',str(O/'code/doublets_full.R'),str(d)],d/'DOUBLETS_COMPLETE')
  run(sm+'_infercnv',[R,'--vanilla',str(O/'code/native_cnv.R'),'infercnv',sm],O/'cnv/infercnv'/sm/'COMPLETE')
 dump('JOB_FINISHED.json',{'exit_code':0,'status':'computation_complete_scientific_review_pending','jobs':jobs,'new_formal_patients':0})
 dump('STATUS.json',{'status':'computation_complete_scientific_review_pending','new_formal_patients':0})
except Exception:
 dump('JOB_FINISHED.json',{'exit_code':1,'status':'failed','traceback':traceback.format_exc(),'jobs':jobs,'new_formal_patients':0})
 dump('STATUS.json',{'status':'failed','new_formal_patients':0});raise
