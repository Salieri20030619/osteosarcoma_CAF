from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,subprocess,signal,time,concurrent.futures,threading,fcntl,traceback,hashlib,re
import numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase95_NCI_expanded_reanalysis_20261016'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
LOCK=threading.Lock();results=[]
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8');t.replace(f)
def stage(sm,name,cmd):
 dump(O/'resources'/(sm+'_STAGE.json'),{'sample':sm,'stage':name,'utc':now()})
 with (O/'logs'/(sm+'_'+name+'.log')).open('w',encoding='utf-8') as log:
  env=os.environ.copy();e=str(Path(R).parent.parent);env.update(PATH=e+'/bin:'+env.get('PATH',''),LD_LIBRARY_PATH=e+'/lib:'+env.get('LD_LIBRARY_PATH',''))
  proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
  try:rc=proc.wait(timeout=16*3600)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=60)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   rc=124
 (O/'logs'/(sm+'_'+name+'.exit_code')).write_text(str(rc)+'\n',encoding='utf-8');return rc
def worker(sm):
 d=O/'inputs'/sm
 while int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))<128*1024*1024:
  dump(O/'resources'/(sm+'_MEMORY_WAIT.json'),{'utc':now(),'required_GiB':128});time.sleep(900)
 assert stage(sm,'prepare',[PY,str(O/'code/prepare.py'),sm])==0 and (d/'PREPARATION_COMPLETE').exists()
 assert stage(sm,'doublets',[R,'--vanilla',str(O/'code/doublets.R'),str(d)])==0 and (d/'DOUBLETS_COMPLETE').exists()
 c=pd.read_csv(d/'cells.tsv',sep='\t',encoding='utf-8');db=pd.read_csv(d/'full_scDblFinder_calls.tsv',sep='\t',encoding='utf-8').set_index('cell_id')
 assert db.index.is_unique and set(db.index)==set(c.cell_id)
 c['full_scDblFinder_class']=c.cell_id.map(db.scDblFinder_class);c['full_scDblFinder_score']=c.cell_id.map(db.scDblFinder_score)
 c['reference_used']=False;c['heldout_clean']=False
 host=c.Ann_Level0.eq('Host') & c.full_scDblFinder_class.eq('singlet')
 masks={'immune':host & c.Ann_Level1.isin(['Immune_Myeloid','Immune_Lymphoid']) & c.provisional_lineage.isin(['T_NK','B_Plasma','Myeloid','Mast']), 'endo':host & c.Ann_Level2.eq('Endothelial_cells') & c.provisional_lineage.eq('Endothelial')}
 reference_plan=[]
 for role,mask in masks.items():
  ids=c.index[mask].to_numpy();rng=np.random.default_rng(int(hashlib.sha256((sm+'|'+role+'|20261016').encode()).hexdigest()[:8],16));ids=rng.permutation(ids)
  # Coarse immune/endo strata are fixed before predictions. Scarce strata remain visible as coverage-limited.
  if len(ids)>=50:
   nh=min(50,len(ids)//5);held=ids[:nh];train=ids[nh:nh+min(200,len(ids)-nh)]
   c.loc[train,['reference','reference_used']]=True;c.loc[held,['heldout','heldout_clean']]=True
   c.loc[train,'reference_role']='training_'+role;c.loc[held,'reference_role']='holdout_'+role
  else:nh=0;train=[]
  reference_plan.append({'stratum':role,'available_clean_author_agreeing':len(ids),'training':len(train),'heldout':nh,'coverage':'assessed' if nh else 'insufficient_pool_not_assessed'})
 c.to_csv(d/'cells.tsv',sep='\t',index=False,encoding='utf-8');dump(d/'REFERENCE_PLAN.json',{'seed_recipe':'SHA256 sample|stratum|20261016;80/20 capped200/50','strata':reference_plan,'references':int(c.reference_used.sum()),'heldout':int(c.heldout_clean.sum())})
 eligible=c.reference_used.sum()>=30 and c.heldout_clean.sum()>=20
 rc=stage(sm,'copykat',['bash',str(O/'code/run_copykat.sh'),sm]) if eligible else None
 marker=(O/'copykat'/sm/'COMPLETE').exists();valid=False;checks=[];actual_refs=0
 if rc==0 and marker:
  calls=pd.read_csv(O/'copykat'/sm/'cell_results.tsv',sep='\t',encoding='utf-8').set_index('cell_id');assert calls.index.is_unique and set(calls.index)==set(c.cell_id)
  c['native_copykat_call']=c.cell_id.map(calls.native_copykat_call)
  log=(O/'logs'/(sm+'_copykat.log')).read_text(encoding='utf-8');nums=re.findall(r'(\d+) known normal cells found in dataset',log)
  assert len(nums)==1;actual_refs=int(nums[0])
  for role,h in c[c.heldout_clean].groupby('reference_role'):
   ev=h.native_copykat_call.isin(['diploid','aneuploid']);n=int(ev.sum());normal=int(h.native_copykat_call.eq('diploid').sum())
   checks.append({'stratum':role,'heldout':len(h),'evaluable':n,'diploid':normal,'fraction':normal/n if n else None,'pass':bool(n>=10 and normal/n>=.9) if n else False})
  valid=actual_refs>=30 and bool(checks) and all(v['pass'] for v in checks)
 else:c['native_copykat_call']='not_assessed_insufficient_references' if not eligible else 'native_execution_failed'
 ar=stage(sm,'ambient',[R,'--vanilla',str(O/'code/ambient.R'),str(d)])
 if ar==0 and (d/'AMBIENT_COMPLETE').exists():
  am=pd.read_csv(d/'ambient_evidence.tsv',sep='\t',encoding='utf-8');assert am.cell_id.is_unique and set(am.cell_id)==set(c.cell_id);c=c.merge(am,on='cell_id',validate='one_to_one')
 c['reference_validation_pass']=valid
 c['mechanical_expanded_candidate']=c.provisional_lineage.eq('Fibroblast') & c.full_scDblFinder_class.eq('singlet') & c.native_copykat_call.eq('diploid') & valid
 c['final_expanded_admission']=False;c['manual_cluster_review_pending']=True
 c.to_csv(O/'tables'/(sm+'_all_cell_evidence.tsv.gz'),sep='\t',encoding='utf-8',index=False)
 c.groupby(['cluster','Ann_Level0','Ann_Level2','provisional_lineage','native_copykat_call','full_scDblFinder_class']).size().rename('cells').reset_index().to_csv(O/'tables'/(sm+'_cluster_evidence.tsv'),sep='\t',encoding='utf-8',index=False)
 result={'sample':sm,'utc':now(),'status':'native_completed' if rc==0 and marker else 'insufficient_references' if not eligible else 'native_failed','copykat_exit':rc,'copykat_marker':marker,'ambient_exit':ar,'cells':len(c),'author_CAF':int(c.author_CAF.sum()),'fibroblast_winner':int(c.provisional_lineage.eq('Fibroblast').sum()),'mechanical_candidates':int(c.mechanical_expanded_candidate.sum()),'reference_pass':bool(valid),'actual_native_references':actual_refs,'reference_checks':checks,'reference_plan':reference_plan,'putative_patient_group':str(c.putative_patient_group.iloc[0]),'independent_patient_verified':False,'manual_review_pending':True}
 dump(O/'resources'/(sm+'_COMPLETE.json'),result);return result
export=json.loads((O/'resources/EXPORT_JOB_FINISHED.json').read_text(encoding='utf-8'));assert export['exit_code']==0 and export['completion_marker']
manifest=json.loads((O/'resources/EXPORTED_INPUTS.json').read_text(encoding='utf-8'));samples=[v['sample'] for v in manifest.values()];assert len(samples)==19
dump(O/'STATUS.json',{'utc':now(),'status':'reanalysis_running','pid':os.getpid(),'samples':19,'manual_review_pending':True})
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 futures={ex.submit(worker,sm):sm for sm in samples}
 for fu in concurrent.futures.as_completed(futures):
  sm=futures[fu]
  try:r=fu.result()
  except Exception:r={'sample':sm,'status':'failed','utc':now(),'traceback':traceback.format_exc()}
  results.append(r);dump(O/'SAMPLE_PROGRESS.json',{'utc':now(),'completed':len(results),'expected':19,'samples':results,'manual_review_pending':True});print('SAMPLE_FINISHED',json.dumps(r),flush=True)
ok=len(results)==19 and all(r['status'] in ['native_completed','insufficient_references'] for r in results)
dump(O/'COMPUTATION_FINISHED.json',{'utc':now(),'exit_code':0 if ok else 1,'samples':results,'manual_review_pending':True,'independent_patient_count':None})
dump(O/'STATUS.json',{'utc':now(),'status':'computation_complete_manual_review_pending' if ok else 'execution_failures_require_review','completed':len(results),'expected':19})
