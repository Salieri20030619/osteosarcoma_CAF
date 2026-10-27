from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,subprocess,signal,time,shutil,concurrent.futures,threading,select,fcntl,traceback
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase92_copykat_validation_20261015'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript';LOCK=threading.Lock();results=[]
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(f,v):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');t.replace(f)
def run(cmd,log):
 with log.open('w',encoding='utf-8') as f:
  p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:return p.wait(timeout=16*3600)
  except subprocess.TimeoutExpired:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=60)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   return 124
def worker(sm):
 src=P/('phase91_lung_patient_validation_20261015' if sm.startswith('GSE') else 'phase90_patient_validation_20261015')/'inputs'/sm
 d=O/'inputs'/sm;d.mkdir(exist_ok=True)
 for name in ['genes.tsv','cells.tsv']:shutil.copy2(src/name,d/name)
 if not (d/'counts.mtx.gz').exists():(d/'counts.mtx.gz').symlink_to(src/'counts.mtx.gz')
 cells=pd.read_csv(d/'cells.tsv',sep='\t',encoding='utf-8')
 source_exit=(P/'phase91_lung_patient_validation_20261015/logs/doublets.exit_code') if sm.startswith('GSE') else (P/'phase90_patient_validation_20261015/logs'/(sm+'_doublets.exit_code'))
 if (src/'DOUBLETS_COMPLETE').exists() and source_exit.exists() and source_exit.read_text(encoding='utf-8').strip()=='0':
  shutil.copy2(src/'full_scDblFinder_calls.tsv',d/'full_scDblFinder_calls.tsv');doublet_origin='existing_completed_full_input_run'
 else:
  rc=run([R,'--vanilla',str(O/'code/doublets_full.R'),str(d)],O/'logs'/(sm+'_doublets.log'))
  (O/'logs'/(sm+'_doublets.exit_code')).write_text(str(rc)+'\n',encoding='utf-8');assert rc==0 and (d/'DOUBLETS_COMPLETE').exists()
  doublet_origin='independent_full_input_run_same_seed_recipe'
 db=pd.read_csv(d/'full_scDblFinder_calls.tsv',sep='\t',encoding='utf-8').set_index('cell_id');assert db.index.is_unique and set(db.index)==set(cells.cell_id)
 cells['full_scDblFinder_class']=cells.cell_id.map(db.scDblFinder_class);assert cells.full_scDblFinder_class.isin(['singlet','doublet']).all()
 cells['reference_used']=cells.reference&cells.full_scDblFinder_class.eq('singlet')
 cells['heldout_clean']=cells.heldout&cells.full_scDblFinder_class.eq('singlet')
 cells['reference_role']=cells.get('reference_role',cells.get('phase88_cnv_role'))
 cells.to_csv(d/'cells.tsv',sep='\t',index=False,encoding='utf-8')
 if cells.reference_used.sum()<30 or cells.heldout_clean.sum()<20:return {'sample':sm,'status':'ineligible_clean_reference','new_formal_patients':0}
 while int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))<128*1024*1024:
  dump(O/'resources'/(sm+'_MEMORY_WAIT.json'),{'utc':now(),'threshold_GiB':128});time.sleep(900)
 dump(O/'resources'/(sm+'_RUNNING.json'),{'utc':now(),'sample':sm,'doublet_origin':doublet_origin,'references':int(cells.reference_used.sum()),'heldout':int(cells.heldout_clean.sum())})
 rc=run(['bash',str(O/'code/run_copykat.sh'),sm],O/'logs'/(sm+'_copykat.log'))
 (O/'logs'/(sm+'_copykat.exit_code')).write_text(str(rc)+'\n',encoding='utf-8')
 complete=(O/'copykat'/sm/'COMPLETE').exists()
 result={'sample':sm,'exit_code':rc,'completion_marker':complete,'doublet_origin':doublet_origin,'utc':now(),'new_formal_patients':0}
 if rc!=0 or not complete:result['status']='failed';return result
 calls=pd.read_csv(O/'copykat'/sm/'cell_results.tsv',sep='\t',encoding='utf-8').set_index('cell_id');assert calls.index.is_unique and set(calls.index)==set(cells.cell_id)
 cells['native_copykat_call']=cells.cell_id.map(calls.native_copykat_call)
 checks=[]
 for label,g in cells[cells.heldout_clean].groupby('reference_role'):
  ev=g.native_copykat_call.isin(['diploid','aneuploid']);n=int(ev.sum());normal=int(g.native_copykat_call.eq('diploid').sum())
  checks.append({'group':label,'cells':len(g),'evaluable':n,'diploid':normal,'diploid_fraction':normal/n if n else None,'pass':bool(n>=10 and normal/n>=.9) if n else False})
 ok=bool(checks) and all(x['pass'] for x in checks)
 cells['copykat_reference_pass']=ok;cells['formal_admission']=False
 cells.to_csv(O/'tables'/(sm+'_copykat_full_evidence.tsv.gz'),sep='\t',index=False,encoding='utf-8')
 result.update(status='complete_reference_pass' if ok else 'complete_reference_failed',reference_checks=checks,reference_pass=ok,native_calls=cells.native_copykat_call.value_counts().to_dict(),full_cells=len(cells))
 dump(O/'resources'/(sm+'_COMPLETE.json'),result);return result
def lane(samples,wait_lung=False):
 if wait_lung:
  pid=json.loads((P/'phase91_lung_patient_validation_20261015/LAUNCH.json').read_text(encoding='utf-8'))['pid']
  try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 for sm in samples:
  try:r=worker(sm)
  except Exception:r={'sample':sm,'status':'failed','traceback':traceback.format_exc(),'utc':now()}
  with LOCK:results.append(r);dump(O/'SAMPLE_PROGRESS.json',{'updated_utc':now(),'samples':results,'scientific_admission_pending':True})
  print('SAMPLE_FINISHED',json.dumps(r),flush=True)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
dump(O/'STATUS.json',{'status':'full_input_native_CopyKAT_running','pid':os.getpid(),'utc':now(),'scientific_admission_pending':True})
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 fs=[ex.submit(lane,['GSE234187__GSM7453693','NCC__BS7','NCC__BS1','NCC__BS12','NCC__BS14']),ex.submit(lane,['NCC__BS9','NCC__BS5','NCC__BS4','NCC__BS3','NCC__BS2'],True)]
 for f in fs:f.result()
ok=len(results)==10 and all(r['status'].startswith('complete_') or r['status'].startswith('ineligible_') for r in results)
dump(O/'JOB_FINISHED.json',{'exit_code':0 if ok else 1,'utc':now(),'results':results,'scientific_admission_pending':True,'new_formal_patients':0})
dump(O/'STATUS.json',{'status':'computation_complete_review_pending' if ok else 'tool_failures_review_required','utc':now(),'new_formal_patients':0})
