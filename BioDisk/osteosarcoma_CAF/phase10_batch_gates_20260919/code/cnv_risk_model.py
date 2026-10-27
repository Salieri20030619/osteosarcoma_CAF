from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import ast,json
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse as sp
from sklearn.metrics import adjusted_rand_score,silhouette_score
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';SEED=20260919
# Reuse the exact model functions without executing the primary pipeline again.
module=ast.parse((O/'code/batch_models.py').read_text());funcs=[n for n in module.body if isinstance(n,ast.FunctionDef) and n.name in ['normalize','pca','harmonize','evaluate']];exec(compile(ast.Module(body=funcs,type_ignores=[]),'<frozen_model_functions>','exec'))
a=ad.read_h5ad(O/'objects/unified_strict_integrated.h5ad');old=ad.read_h5ad(P/'phase9_joint_20260919/objects/joint_CAF_integrated.h5ad');extra=pd.read_csv(O/'tables/CNV_risk_0.3_additional_candidates.tsv.gz',sep='\t',index_col='cell_id');assert len(extra)==91 and extra.old_joint_member.all();assert set(extra.index).isdisjoint(a.obs_names)
ix=old.obs_names.get_indexer(extra.index);vi=old.var_names.get_indexer(a.var_names);assert (ix>=0).all() and (vi>=0).all();raw=sp.vstack([a.layers['counts'],old.layers['counts'][ix][:,vi]]).tocsr();names=a.obs_names.append(extra.index)
e=pd.read_csv(O/'tables/all_cell_unified_evidence.tsv.gz',sep='\t',index_col='cell_id',low_memory=False);ob=e.loc[names].copy();ob.index=pd.Index(ob.index.to_numpy(dtype=object),dtype=object)
for c in ob.select_dtypes(include=['object','string']):ob[c]=ob[c].fillna('').astype(object)
b=ad.AnnData(normalize(raw,ob.full_library_counts),obs=ob,var=a.var.copy());b.layers['counts']=raw;b.obs['CNV_risk_0_3_eligible']=True
assert b.obs.dual_low.all() and b.obs.ambient_checked.all() and ~b.obs.independent_doublet_excluded.any();assert (b.obs.cluster_any_high_fraction<.3).all()
metrics=[];support=[];overlap=[];z,_=pca(b);evaluate(b,harmonize(z,b),'CNV_risk_0_3',umap=True);rr=a.obs.unified_cluster.astype(str);ll=b.obs.loc[a.obs_names,'CNV_risk_0_3'].astype(str);metrics[0]['ARI_vs_unified_primary']=adjusted_rand_score(rr,ll)
for c in sorted(rr.unique(),key=int):
 u=rr==c;best=max((int((u&(ll==v)).sum())/int((u|(ll==v)).sum()),v) for v in ll.unique());overlap.append({'analysis':'CNV_risk_0_3','primary_cluster':c,'primary_cells_shared':int(u.sum()),'best_cluster':best[1],'Jaccard':best[0]})
b.write_h5ad(O/'objects/CNV_risk_0_3_integrated.h5ad',compression='gzip')
for name,rows in [('batch_model_comparison.tsv',metrics),('cluster_support.tsv',support),('model_cluster_overlap.tsv',overlap)]:
 d=pd.read_csv(O/'tables'/name,sep='\t');d=d[d.analysis!='CNV_risk_0_3'];pd.concat([d,pd.DataFrame(rows)],ignore_index=True).to_csv(O/'tables'/name,sep='\t',index=False)
s=json.loads((O/'STATUS.json').read_text());s['batch_models']=len(pd.read_csv(O/'tables/batch_model_comparison.tsv',sep='\t'));s['CNV_risk_0_3_sensitivity_cells']=b.n_obs;(O/'STATUS.json').write_text(json.dumps(s,indent=2));(O/'CNV_RISK_SENSITIVITY_COMPLETE').write_text('CNV_RISK_SENSITIVITY_COMPLETE\n');print(json.dumps(metrics[0]),flush=True)