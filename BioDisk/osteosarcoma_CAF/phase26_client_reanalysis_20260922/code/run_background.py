from common import *
import threading,fcntl,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed

def native_job(method,sm):
 d=O/'cnv'/method/sm;d.mkdir(parents=True,exist_ok=True);log=O/'logs'/f'{method}_{sm}.log'
 if (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists():return {'method':method,'status':'complete_cached','exit_code':0,'path':str(d)}
 with log.open('a') as out:rc=subprocess.run(['prlimit','--as=103079215104',R,'--vanilla',str(O/'code/native_cnv.R'),method,sm],stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT).returncode
 (d/'exit_code').write_text(str(rc)+'\n');ok=rc==0 and (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists()
 return {'method':method,'status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(d),'failure':(d/'FAILED').read_text().strip() if (d/'FAILED').exists() else ''}

def combine(row,jobs):
 sm=row['sample'];d=O/'inputs'/sm;cells=read(d/'cells.tsv');ids=pd.Index(cells.cell_id);ev=cells[['cell_id','n_counts','reference']].copy();ev['sample']=sm;ev['dataset']=row['dataset'];ev['copykat_call']='not_evaluated';ev['infercnv_call']='not_evaluated';ev['infercnv_source']='not_evaluated'
 for j in jobs:
  if not j['status'].startswith('complete'):continue
  tab=read(Path(j['path'])/'cell_results.tsv').set_index('cell_id');ev[j['method']+'_call']=tab['call'].reindex(ids).fillna('not_returned_by_software').to_numpy()
  if j['method']=='infercnv':ev['infercnv_source']='native_i6_HMM'
 prov=read_json_local(d/'provider_evidence.json')
 if prov.get('infercnv_status')=='success' and (d/'provider_infercnv_table.tsv.gz').exists():
  t=read(d/'provider_infercnv_table.tsv.gz');cols=[c for c in t if c.startswith('has_cnv_chr')];assert cols and 'barcodes' in t
  t.index=sm+'_'+t.barcodes.astype(str).str.strip();assert t.index.is_unique
  vals=t[cols].apply(pd.to_numeric,errors='coerce');valid=vals.notna().all(axis=1);calls=pd.Series('not_evaluated',index=t.index);calls.loc[valid]=np.where(vals.loc[valid].eq(1).any(axis=1),'alteration_called','no_alteration_called');ev['infercnv_call']=calls.reindex(ids).fillna('not_evaluated').to_numpy();ev['infercnv_source']='official_provider_HMM_calls'
 ck=ev.copykat_call;ic=ev.infercnv_call;lowck=ck.eq('diploid');highck=ck.eq('aneuploid');lowic=ic.eq('no_alteration_called');highic=ic.eq('alteration_called')
 ev['CNV_evidence']='unassessed';ev.loc[lowck|lowic,'CNV_evidence']='single_method_no_alteration';ev.loc[highck|highic,'CNV_evidence']='single_method_alteration';ev.loc[lowck&lowic,'CNV_evidence']='concordant_no_alteration';ev.loc[highck&highic,'CNV_evidence']='concordant_alteration';ev.loc[(lowck&highic)|(highck&lowic),'CNV_evidence']='method_conflict'
 a=ad.read_h5ad(d/'full_called_counts.h5ad',backed='r')
 for c in ['provider_consensus_celltype_annotation','historical_broad_lineage','provider_is_infercnv_reference']:
  if c in a.obs:ev[c]=a.obs[c].astype(str).reindex(ids).to_numpy()
 a.file.close();ev['CAF_identity_review_status']='CNV_EVIDENCE_ONLY_NOT_CONFIRMED_CAF'
 ev.to_csv(O/'tables'/f'{sm}_identity_evidence.tsv.gz',sep='\t',index=False)
 return ev.CNV_evidence.value_counts().to_dict()

def read_json_local(p):return json.loads(p.read_text()) if p.exists() else {}
def run_sample(row,qa):
 sm=row['sample'];d=O/'inputs'/sm;write_json(d/'JOB_STATE.json',{'sample':sm,'status':'CNV_RUNNING','started':now()});jobs=[native_job('copykat',sm)]
 reused=qa['provider_infercnv_status']=='success' and (d/'provider_infercnv_table.tsv.gz').exists()
 if not reused:jobs.append(native_job('infercnv',sm))
 counts=combine(row,jobs);status='CNV_EVIDENCE_READY' if all(j['status'].startswith('complete') for j in jobs) else 'CNV_EVIDENCE_WITH_SOFTWARE_FAILURE'
 result={'sample':sm,'dataset':row['dataset'],'status':status,'input_cells':qa['cells'],'official_infercnv_reused':reused,'jobs':jobs,'evidence_counts':counts,'finished':now()};write_json(d/'JOB_STATE.json',result);return result

lock=(O/'runner.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
manifest=read(O/'tables/analysis_source_manifest.tsv').fillna('');rows=manifest.to_dict('records');assert len(rows)==107
# Interleave cohorts for scheduling only; all records remain scheduled and participate in final accounting.
groups=[g.to_dict('records') for _,g in manifest.groupby('dataset',sort=True)];ordered=[]
while any(groups):
 for group in groups:
  if group:ordered.append(group.pop(0))
state={'status':'RUNNING_INPUTS_AND_CNV','started_at':now(),'source_records':107,'prepared':0,'finished':0,'input_failures':0,'concurrency':2,'current_input':None,'formal_CAF_clustering_started':False};guard=threading.RLock();results={};futures=[]
def snapshot():
 with guard:
  state['updated_at']=now();state['finished']=len(results);write_json(O/'STATUS.json',state);write_json(O/'tables/job_outcomes.json',list(results.values()))
def callback(f):
 sm=f.sample
 try:r=f.result()
 except Exception as e:r={'sample':sm,'status':'PROCESSING_FAILED','error':str(e),'finished':now()}
 with guard:results[sm]=r;snapshot()
rc=1
try:
 (O/'RUN_STARTED').write_text(str(os.getpid())+'\n');snapshot()
 with ThreadPoolExecutor(max_workers=2) as pool:
  for row in ordered:
   sm=row['sample']
   with guard:state['current_input']=sm;snapshot()
   try:
    qa=prepare(row)
    with guard:state['prepared']+=1;snapshot()
    f=pool.submit(run_sample,row,qa);f.sample=sm;f.add_done_callback(callback);futures.append(f);print('INPUT_READY',sm,qa['cells'],flush=True)
   except Exception as e:
    with guard:state['input_failures']+=1;results[sm]={'sample':sm,'dataset':row['dataset'],'status':'INPUT_FAILED','error':str(e)};snapshot()
    print('INPUT_FAILED',sm,type(e).__name__,str(e),flush=True)
   gc.collect()
  with guard:state['current_input']=None;state['status']='WAITING_FOR_NATIVE_CNV';snapshot()
  for f in as_completed(futures):
   try:f.result()
   except Exception:pass
 assert len(results)==107
 summaries=[]
 for row in rows:
  r=results[row['sample']];summaries.append({'sample':row['sample'],'dataset':row['dataset'],'status':r['status'],'input_cells':r.get('input_cells'),**r.get('evidence_counts',{})})
 pd.DataFrame(summaries).to_csv(O/'tables/all107_sample_outcomes.tsv',sep='\t',index=False)
 bad=[r for r in results.values() if r['status'] in ['INPUT_FAILED','PROCESSING_FAILED']]
 qa={'finished_at':now(),'scheduled_sample_records':107,'recorded_outcomes':107,'input_or_processing_failures':len(bad),'software_failure_samples':sum(r['status']=='CNV_EVIDENCE_WITH_SOFTWARE_FAILURE' for r in results.values()),'formal_CAF_clustering_started':False,'next_action':'Review all-sample CNV and lineage identity, fix unresolved technical failures, then proceed to CAF clustering and myCAF reclustering without invented thresholds','all_software_success_implied':False}
 write_json(O/'delivery/CNV_STAGE_QA.json',qa);(O/'CNV_STAGE_FINISHED').write_text('ALL_107_SOURCE_RECORDS_HAVE_OUTCOMES\n');state['status']='INPUT_REPAIR_REQUIRED' if bad else 'CNV_IDENTITY_REVIEW_REQUIRED';rc=1 if bad else 0;snapshot()
except Exception as e:
 state['status']='RUNNER_FAILED';state['error']=str(e);snapshot();traceback.print_exc()
finally:
 (O/'exit_code').write_text(str(rc)+'\n');(O/'RUN_FINISHED').write_text(now()+'\n')
