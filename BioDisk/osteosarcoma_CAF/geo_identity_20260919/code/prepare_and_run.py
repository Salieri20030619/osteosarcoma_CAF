from pathlib import Path
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='2'
import json,datetime,gzip,subprocess,concurrent.futures,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'geo_identity_20260919'; Q=P/'geo_preprocess_20260918'
for d in ['inputs','decontx','logs','tables','figures','objects']: (O/d).mkdir(parents=True,exist_ok=True)
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
manifest=pd.read_csv(Q/'manifest.tsv',sep='\t').fillna('')
manifest.to_csv(O/'tables/sample_provenance.tsv',sep='\t',index=False)
rules={'version':'GEO_identity_20260919_v1','frozen_utc':stamp(),'base':'phase7_review_states_20260917/code/integrate_identity.py',
 'strict':'FibLike, >=2 Fib genes, n_genes>=500, both CNV low, each used reference/holdout stratum >=10 evaluable and >=90% normal, no CTA/competing lineage, cluster concordant-high fraction <0.2',
 'additional_identity_guard':'exclude clusters with >=20% any-method elevated among jointly binary-evaluable cells; require no strong Osteo program (>=3 genes AND >0.2% counts)',
 'ambient':'per-sample decontX, broad reviewed lineage labels, default delta/estimateDelta, no empty droplet background; counts>=1 main corrected detection, >0 sensitivity; original evidence not rescued by correction',
 'corrected_support':'raw strict + same Fib>=2 and competing/CTA gates on corrected counts; not a new CNV analysis',
 'FLEX':'operational candidates only; separate technology flag, no independently validated CNV sensitivity',
 'patient_policy':'no pooled cross-cohort N or clinical inference; retain provenance and overlap-unresolved flags',
 'no_biological_truth_claim':True,'no_inferential_tests':True,'original_labels_unchanged':True}
(O/'decision_rules.json').write_text(json.dumps(rules,indent=2))
jobs=[]; hashes=[]
for row in manifest.to_dict('records'):
 sm=row['sample']; evidence=pd.read_csv(P/'geo_cnv_20260918/review'/(sm+'_cells.tsv.gz'),sep='\t',dtype={'cluster':str})
 if not evidence.lineage_phase3.eq('FibLike').any(): continue
 src=Q/'accepted'/sm/'analysis_ready.h5ad'; a=ad.read_h5ad(src)
 assert a.obs_names.tolist()==evidence.cell_id.tolist()
 d=O/'inputs'/sm;d.mkdir(exist_ok=True)
 x=a.layers['counts'].tocsr(); assert (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
 with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
 pd.DataFrame({'cell_id':a.obs_names,'cluster':evidence.lineage_phase3.to_numpy()}).to_csv(d/'cells.tsv',sep='\t',index=False)
 pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.to_numpy()}).to_csv(d/'genes.tsv',sep='\t',index=False)
 with open(src,'rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
 hashes.append({'path':str(src),'sha256':digest})
 (d/'COMPLETE').write_text('COUNTS_AND_CELL_ORDER_VERIFIED\n');jobs.append(sm)
 print('EXPORTED',sm,a.shape,flush=True)
 del a,x
pd.DataFrame(hashes).to_csv(O/'tables/accepted_input_sha256.tsv',sep='\t',index=False)
(O/'STATUS.json').write_text(json.dumps({'status':'AMBIENT_SENSITIVITY_RUNNING','started_utc':stamp(),'samples':jobs},indent=2))
def run(sm):
 with open(O/'logs'/(sm+'.log'),'w') as f:
  rc=subprocess.call(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/run_decontx.R'),sm],stdout=f,stderr=subprocess.STDOUT)
 (O/'logs'/(sm+'.exit_code')).write_text(str(rc)+'\n')
 print('DECONTX_EXIT',sm,rc,flush=True)
 return {'sample':sm,'exit_code':rc,'complete':(O/'decontx'/sm/'COMPLETE').exists()}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex: results=list(ex.map(run,jobs))
pd.DataFrame(results).to_csv(O/'tables/decontx_jobs.tsv',sep='\t',index=False)
assert all(x['exit_code']==0 and x['complete'] for x in results)
(O/'AMBIENT_COMPLETE').write_text('ALL_SAMPLE_CORRECTIONS_COMPLETE\n')
(O/'STATUS.json').write_text(json.dumps({'status':'AMBIENT_COMPLETE_IDENTITY_REVIEW_PENDING','completed_utc':stamp(),'samples':len(jobs)},indent=2))
print('AMBIENT_COMPLETE',len(jobs),flush=True)
