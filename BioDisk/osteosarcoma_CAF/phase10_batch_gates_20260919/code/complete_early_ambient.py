from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,anndata as ad,numpy as np,json,gzip,subprocess,concurrent.futures,hashlib
from scipy.io import mmwrite
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919'
e=pd.read_csv(P/'phase3_identity_20260914/identity/all_cell_identity_evidence.tsv.gz',sep='\t',low_memory=False);jobs=[]
for sm,g in e.groupby('sample',sort=True):
 a=ad.read_h5ad(P/'OS_CAF_phase2/results'/sm/'annotated.h5ad');g=g.set_index('cell_id').loc[a.obs_names];x=a.layers['counts'].tocsr();assert np.array_equal(np.asarray(x.sum(1)).ravel(),g.n_counts.to_numpy())
 d=O/'inputs'/sm;d.mkdir(exist_ok=True)
 if not (d/'COMPLETE').exists():
  with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
  pd.DataFrame({'cell_id':a.obs_names,'cluster':g.lineage_phase3.to_numpy()}).to_csv(d/'cells.tsv',sep='\t',index=False);pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.to_numpy()}).to_csv(d/'genes.tsv',sep='\t',index=False);(d/'COMPLETE').write_text('EXACT_FULL_COUNTS_VERIFIED\n')
 jobs.append(sm);print('EXPORTED',sm,a.shape,flush=True)
script=(P/'geo_identity_20260919/code/run_decontx.R').read_text().replace('/geo_identity_20260919','/phase10_batch_gates_20260919');(O/'code/run_decontx.R').write_text(script)
def run(sm):
 if (O/'decontx'/sm/'COMPLETE').exists() and (O/'logs'/(sm+'.exit_code')).exists() and (O/'logs'/(sm+'.exit_code')).read_text().strip()=='0':return {'sample':sm,'exit_code':0,'complete':True,'reused':True}
 with open(O/'logs'/(sm+'.log'),'w') as f:rc=subprocess.call(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/run_decontx.R'),sm],stdout=f,stderr=subprocess.STDOUT)
 (O/'logs'/(sm+'.exit_code')).write_text(str(rc)+'\n');print('DECONTX_EXIT',sm,rc,flush=True);return {'sample':sm,'exit_code':rc,'complete':(O/'decontx'/sm/'COMPLETE').exists()}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:res=list(pool.map(run,jobs))
pd.DataFrame(res).to_csv(O/'tables/early_decontx_jobs.tsv',sep='\t',index=False);assert all(r['exit_code']==0 and r['complete'] for r in res)
(O/'EARLY_AMBIENT_COMPLETE').write_text('EARLY_SEVEN_SAMPLE_CORRECTIONS_COMPLETE\n');print('EARLY_AMBIENT_COMPLETE',flush=True)