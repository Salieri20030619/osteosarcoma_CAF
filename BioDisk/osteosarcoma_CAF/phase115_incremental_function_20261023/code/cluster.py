from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='2'
import json,re
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023'
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
X=sparse.load_npz(O/'objects/log1p_cp10k_shared.npz');m=pd.read_csv(O/'tables/cell_mapping.tsv.gz',sep='\t',low_memory=False);genes=pd.read_csv(O/'tables/shared_genes.tsv',sep='\t').symbol.to_numpy();ex=set(json.loads((O/'resources/MAPPING_EXCLUDED_GENES.json').read_text()))
allowed=np.array([g not in ex and not re.match('^(MT-|RPL|RPS)',g) and g not in {'XIST','TSIX','MALAT1','DDX3Y','KDM5D','UTY','EIF1AY','RPS4Y1','RPS4Y2','ZFY','USP9Y'} for g in genes])
dump('resources/FUNCTION_CLUSTER_PROTOCOL.json',dict(reason='Heldout transfer diagnostics do not support reliable F2/F2C3 identity transfer; complementary within-source exploratory grouping without tested functional genes',original_mapping_preserved=True,query_sources=['NCC','NCI'],fit_cells='Only singlets; doublets retained in registry and assigned nearest singlet cluster by PCA centroid for descriptive all-cell sensitivity',representation='Existing full-library log1pCP10k; source-specific seurat HVG2000 excluding all functional/lineage/target gene endpoints and technical genes; scale max10; PCA30',neighbors=15,primary_leiden_resolution=.5,primary_seed=20261023,sensitivity_resolution=.8,sensitivity_seed=20261024,no_UMAP_or_figures=True,tests='Same donor/site within-source comparisons, original programs; class/state/cluster FDR families separate',label_policy='NCC_Q0 etc are new computational groups, not renamed F classes or established CAF subtypes',novelty='Not explained by old label transfer does not establish a new biological subtype'))
out=[];qa=[];stability=[]
for src in ['NCC','NCI']:
 rows=np.flatnonzero(m.source.eq(src));mr=m.iloc[rows].copy();sing=~mr.doublet.to_numpy(bool);ix=np.flatnonzero(sing)
 a=ad.AnnData(X[rows[sing]][:,allowed].copy(),obs=pd.DataFrame(index=mr.loc[~mr.doublet,'cell_id'].to_numpy()),var=pd.DataFrame(index=genes[allowed]))
 sc.pp.highly_variable_genes(a,n_top_genes=min(2000,a.n_vars),flavor='seurat');sel=a.var.highly_variable.to_numpy();selected=genes[allowed][sel];assert not set(selected)&ex
 raw=X[rows][:,allowed][:,sel].toarray();mean=raw[sing].mean(0);sd=np.maximum(raw[sing].std(0,ddof=1),1e-8);z=(raw-mean)/sd;z=np.minimum(z,10)
 from sklearn.decomposition import PCA
 pca=PCA(n_components=min(30,len(ix)-1),svd_solver='arpack',random_state=20261023);emb=pca.fit_transform(z[sing]);all_emb=pca.transform(z)
 a=ad.AnnData(z[sing]);a.obsm['X_pca']=emb;sc.pp.neighbors(a,n_neighbors=15,n_pcs=emb.shape[1],use_rep='X_pca',random_state=20261023)
 for name,res,seed in [('primary',.5,20261023),('seed_repeat',.5,20261024),('resolution08',.8,20261023)]:sc.tl.leiden(a,resolution=res,random_state=seed,key_added=name)
 labels=a.obs.primary.astype(str).to_numpy();classes=sorted(set(labels));cent=np.array([emb[labels==c].mean(0) for c in classes]);dist=((all_emb[:,None,:]-cent[None,:,:])**2).sum(2);pred=np.array(classes)[dist.argmin(1)];pred[sing]=labels
 mr['cluster_compatible']=[src+'_Q'+v for v in pred];mr['cluster_fit_singlet']=sing;out.append(mr[['cell_id','cluster_compatible','cluster_fit_singlet']])
 for label in classes:
  mask=labels==label;cellrows=mr.loc[~mr.doublet].iloc[np.flatnonzero(mask)];dc=cellrows.donor.value_counts();jr=[]
  for key in ['seed_repeat','resolution08']:
   other=a.obs[key].astype(str).to_numpy();best=max([np.sum(mask&(other==t))/np.sum(mask|(other==t)) for t in set(other)]);jr.append(best)
  stability.append(dict(source=src,cluster=src+'_Q'+label,singlet_cells=int(mask.sum()),source_patients=len(dc),dominant_patient=dc.index[0],dominant_patient_fraction=float(dc.iloc[0]/mask.sum()),patients_with_at_least5=int((dc>=5).sum()),seed_best_jaccard=jr[0],resolution08_best_jaccard=jr[1]))
 save(pd.DataFrame({'symbol':selected}),src+'_cluster_features.tsv');qa.append(dict(source=src,singlets=int(sing.sum()),doublets=int((~sing).sum()),clusters=len(classes),seed_ARI=adjusted_rand_score(labels,a.obs.seed_repeat),resolution08_ARI=adjusted_rand_score(labels,a.obs.resolution08),test_genes_disjoint=True))
 print('CLUSTER',src,qa[-1],flush=True)
save(pd.concat(out),'query_cluster_membership.tsv.gz');save(pd.DataFrame(stability),'query_cluster_stability.tsv');dump('CLUSTER_COMPLETE.json',dict(status='PASS_EXPLORATORY',sources=qa))
