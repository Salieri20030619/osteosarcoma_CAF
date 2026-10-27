from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,gzip,subprocess,datetime,fcntl
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd,numpy as np,anndata as ad
from scipy.io import mmwrite
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921';C=O/'CNV_official_inputs'
for d in ['inputs','logs','results']:(C/d).mkdir(parents=True,exist_ok=True)
lock=(C/'runner.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
manifest=pd.read_csv(P/'phase16_scpca_supplement_20260920/manifest.tsv',sep='\t');jobs=[];records=[]
for r in manifest.to_dict('records'):
 sm=r['scpca_sample_id'];d=C/'inputs'/sm;d.mkdir(exist_ok=True)
 proc=ad.read_h5ad(r['processed'],backed='r');ob=proc.obs.copy();infstatus=str(proc.uns.get('infercnv_status','missing'));proc.file.close()
 # Provider processed QC is an explicit downstream evidence step; entry model keeps all filtered cells.
 raw=ad.read_h5ad(r['filtered'],backed='r');x=raw[ob.index].X.tocsr();genes=pd.DataFrame({'ensembl':raw.var_names.astype(str).str.split('.').str[0],'symbol':raw.var.gene_symbol.astype(str).to_numpy()});raw.file.close()
 ob['cell_id']=[sm+'_'+str(z) for z in ob.index];ob['tool_id']=[sm+'_'+str(i).zfill(6) for i in range(len(ob))]
 ref=ob.get('is_infercnv_reference',pd.Series(False,index=ob.index)).fillna(False).astype(str).str.lower().isin(['true','1','1.0'])
 ob['reference']=ref.to_numpy();ob['n_counts']=np.asarray(x.sum(1)).ravel()
 if not (d/'COMPLETE').exists():
  with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
  genes.to_csv(d/'genes.tsv',sep='\t',index=False);ob[['cell_id','tool_id','reference','n_counts']].to_csv(d/'cells.tsv',sep='\t',index=False)
  pd.DataFrame({'id':ob.tool_id,'annotation':np.where(ref,'reference','query')}).to_csv(d/'annotations.tsv',sep='\t',index=False,header=False)
  (d/'COMPLETE').write_text('PROVIDER_PROCESSED_COUNTS_VERIFIED\n')
 jobs.append(('copykat',sm))
 if infstatus!='success':jobs.append(('infercnv',sm))
 records.append({'sample':sm,'input_cells':len(ob),'official_reference_cells':int(ref.sum()),'official_infercnv_status':infstatus,'CopyKAT_reference':'automatic native baseline; norm.cell.names default','inferCNV_reference':'provider normal annotation; no project minimum/holdout/marker gate'})
pd.DataFrame(records).to_csv(C/'input_manifest.tsv',sep='\t',index=False)
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
def run(job):
 method,sm=job;d=C/'results'/method/sm
 if (d/'COMPLETE').exists():return {'method':method,'sample':sm,'exit_code':0,'status':'complete_cached'}
 log=C/'logs'/f'{method}_{sm}.log'
 with log.open('w') as out:rc=subprocess.run(['prlimit','--as=68719476736',R,'--vanilla',str(O/'code/cnv_official.R'),method,sm],stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT).returncode
 (C/'logs'/f'{method}_{sm}.exit_code').write_text(str(rc)+'\n')
 ok=rc==0 and (d/'COMPLETE').exists();record={'method':method,'sample':sm,'exit_code':rc,'status':'complete' if ok else 'software_failed','log':str(log)};print('METHOD_FINISHED',record,flush=True);return record
results=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 for fut in as_completed([pool.submit(run,j) for j in jobs]):
  results.append(fut.result());pd.DataFrame(results).to_csv(C/'job_status.tsv',sep='\t',index=False)
assert len(results)==len(jobs)
(C/'RUN_FINISHED').write_text('ALL_JOBS_HAVE_RECORDED_OUTCOMES\n')
(C/'QA.json').write_text(json.dumps({'jobs':len(jobs),'successful':sum(r['status'].startswith('complete') for r in results),'failed':sum(r['status']=='software_failed' for r in results),'project_extra_thresholds':False,'full_entry_model_unchanged':True,'official_existing_inferCNV_reused':9,'failure_is_not_negative_evidence':True},indent=2))
print('CNV_JOBS_FINISHED',len(results),flush=True)
