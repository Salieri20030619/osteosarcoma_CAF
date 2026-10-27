from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,scanpy as sc,harmonypy as hm
from sklearn.metrics import adjusted_rand_score,silhouette_score
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();out=p/'integration_review';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'analysis/GSE152048_strict_CAF_states.h5ad');original=a.obs.caf_cluster.astype(str).copy();base=a.obsm['X_pca'].copy();patients=a.obs.patient_id.astype(str).to_numpy()
programs=[c for c in a.obs if c.startswith('program_')]
rows=[];support=[]
def assess(rep,labels,name):
 labels=np.asarray(labels).astype(str)
 neighbors=NearestNeighbors(n_neighbors=16).fit(rep).kneighbors(rep,return_distance=False)[:,1:]
 rows.append({'analysis':name,'clusters':len(set(labels)),'ARI_vs_unintegrated':adjusted_rand_score(original,labels),'patient_silhouette':silhouette_score(rep,patients),'mean_same_patient_neighbor_fraction':float((patients[neighbors]==patients[:,None]).mean())})
 counts=pd.crosstab(a.obs.patient_id,labels);counts.to_csv(out/(name+'_patient_counts.tsv'),sep='\t')
 eligible=((counts>=20)&(counts.sum(axis=1).to_numpy()[:,None]-counts>=20)).sum(axis=0)
 for c in counts:
  z=a.obs.iloc[np.where(labels==c)[0]];v=counts[c]
  row={'analysis':name,'state':c,'cells':int(v.sum()),'patients_ge10':int((v>=10).sum()),'max_patient_fraction':float(v.max()/v.sum()),'eligible_paired_patients_ge20':int(eligible[c])}
  row.update({k:float(z[k].median()) for k in programs});support.append(row)
assess(base,original,'unintegrated')
for theta in [1,2]:
 h=hm.run_harmony(base,a.obs,'patient_id',theta=theta,random_state=20260917,ncores=4,max_iter_harmony=20,verbose=True)
 corrected=np.asarray(h.Z_corr)
 if corrected.shape!=base.shape:corrected=corrected.T
 assert corrected.shape==base.shape and np.isfinite(corrected).all()
 key=f'X_harmony_theta{theta}';a.obsm[key]=corrected
 sc.pp.neighbors(a,use_rep=key,n_neighbors=15,random_state=20260917)
 for resolution in [.3,.5,.8]:
  name=f'harmony_theta{theta}_r{resolution}'
  sc.tl.leiden(a,resolution=resolution,random_state=20260917,key_added=name,flavor='igraph',n_iterations=2,directed=False)
  assess(corrected,a.obs[name],name)
 sc.tl.umap(a,random_state=20260917);a.obsm[f'X_umap_harmony_theta{theta}']=a.obsm['X_umap'].copy()
 fig,axes=plt.subplots(1,3,figsize=(17,5))
 for ax,keylabel in zip(axes,['patient_id',f'harmony_theta{theta}_r0.5','caf_cluster']):
  sc.pl.embedding(a,basis=f'umap_harmony_theta{theta}',color=keylabel,ax=ax,show=False,title=keylabel)
 fig.tight_layout();fig.savefig(out/f'harmony_theta{theta}_comparison.png',dpi=160,bbox_inches='tight');plt.close(fig)
print('COUNTS_UNCHANGED',a.shape)
pd.DataFrame(rows).to_csv(out/'integration_comparison.tsv',sep='\t',index=False)
pd.DataFrame(support).to_csv(out/'state_support_programs.tsv',sep='\t',index=False)
a.write_h5ad(out/'strict_CAF_harmony_review.h5ad',compression='gzip')
(out/'parameters.json').write_text(json.dumps({'method':'Harmony','version':hm.__version__,'batch':'patient_id','theta':[1,2],'resolutions':[.3,.5,.8],'seed':20260917,'original_identity_thresholds_changed':False,'original_counts_changed':False,'interpretation':'batch sensitivity; increased patient mixing alone does not validate a biological subtype; one lung donor confounds patient and site'},indent=2))
(out/'COMPLETE').write_text('HARMONY_REVIEW_COMPUTE_COMPLETE\n');print(pd.DataFrame(rows).to_string(index=False));print('HARMONY_REVIEW_COMPUTE_COMPLETE',flush=True)
