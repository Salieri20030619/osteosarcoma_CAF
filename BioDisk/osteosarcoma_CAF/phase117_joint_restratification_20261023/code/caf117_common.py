from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,re,sys,gc
import numpy as np,pandas as pd
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase117_joint_restratification_20261023';B=R/'phase115_incremental_function_20261023'
for d in ['code','resources','tables','objects','logs','receipts','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(n,v):
 p=O/n;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8');t.replace(p)
def event(stage,status,**kw):
 rec=dict(utc=now(),stage=stage,status=status,**kw)
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(rec,ensure_ascii=False,default=str)+'\n')
 print(json.dumps(rec,ensure_ascii=False,default=str),flush=True);dump('STATUS.json',rec)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def load():return sparse.load_npz(B/'objects/log1p_cp10k_shared.npz'),read('frozen_registry.tsv.gz'),pd.read_csv(B/'tables/shared_genes.tsv',sep='\t').symbol.to_numpy()
def allowed_genes(genes,exclude=()):
 bad=set(exclude)|{'XIST','TSIX','MALAT1','DDX3Y','KDM5D','UTY','EIF1AY','RPS4Y1','RPS4Y2','ZFY','USP9Y'}
 return np.array([not re.match(r'^(MT-|RPL|RPS)',g) and g not in bad for g in genes])
def representation(X,m,genes,fitrows,prefix,exclude=(),fixed_features=None):
 import scanpy as sc,anndata as ad
 from sklearn.decomposition import PCA
 allowed=allowed_genes(genes,exclude);ix=np.flatnonzero(allowed)
 if fixed_features is None:
  a=ad.AnnData(X[fitrows][:,ix].copy(),obs=m.iloc[fitrows][['source']].copy(),var=pd.DataFrame(index=genes[ix]))
  sc.pp.highly_variable_genes(a,n_top_genes=min(3000,len(ix)),batch_key='source',flavor='seurat');chosen=ix[a.var.highly_variable.to_numpy()];del a
 else:
  chosen=pd.Index(genes).get_indexer(fixed_features);assert (chosen>=0).all();assert allowed[chosen].all()
 raw=X[:,chosen].toarray().astype(np.float32);mu=raw[fitrows].mean(0);sd=np.maximum(raw[fitrows].std(0,ddof=1),1e-8);z=np.minimum((raw-mu)/sd,10)
 pca=PCA(n_components=min(50,len(fitrows)-1),svd_solver='arpack',random_state=20261023);pca.fit(z[fitrows]);pc=pca.transform(z).astype(np.float32)
 np.savez_compressed(O/'objects'/(prefix+'_PCA_model.npz'),feature_indices=chosen,features=genes[chosen],mean=mu,sd=sd,components=pca.components_,pca_mean=pca.mean_,explained_variance=pca.explained_variance_ratio_)
 save(pd.DataFrame({'symbol':genes[chosen]}),prefix+'_HVG.tsv');del raw,z;gc.collect();return pc,genes[chosen]
def harmony(pc,meta,theta):
 import harmonypy as hm
 h=hm.run_harmony(pc,meta[['source']].reset_index(drop=True),['source'],theta=theta,max_iter_harmony=20,random_state=20261023,ncores=4,verbose=False)
 z=np.asarray(h.Z_corr);z=z.T if z.shape==pc.T.shape else z;assert z.shape==pc.shape and np.isfinite(z).all();return z.astype(np.float32)
def clusters(emb,res=.5,seed=20261023):
 import scanpy as sc,anndata as ad
 a=ad.AnnData(np.zeros((len(emb),0),dtype=np.float32));a.obsm['rep']=emb
 sc.pp.neighbors(a,n_neighbors=min(15,len(emb)-1),use_rep='rep',random_state=20261023)
 sc.tl.leiden(a,resolution=res,random_state=seed,key_added='labels',flavor='leidenalg',n_iterations=-1)
 return a.obs.labels.astype(str).to_numpy(),a.obsp['connectivities'].tocsr()
def relabel(labels):return np.array(['J'+str(int(v)).zfill(2) for v in labels])
def neighbors(x,k=15):
 from sklearn.neighbors import NearestNeighbors
 return NearestNeighbors(n_neighbors=k+1,n_jobs=4).fit(x).kneighbors(x,return_distance=False)[:,1:]
def project(pc_train,corrected_train,pc_query,train_sources,query_sources):
 from sklearn.neighbors import NearestNeighbors
 pred=pc_query.copy()
 for src in np.unique(query_sources):
  tr=np.flatnonzero(train_sources==src);qr=np.flatnonzero(query_sources==src)
  if not len(tr):continue
  nn=NearestNeighbors(n_neighbors=min(30,len(tr)),n_jobs=4).fit(pc_train[tr]);dist,idx=nn.kneighbors(pc_query[qr]);w=1/np.maximum(dist,1e-4);w=w/w.sum(1,keepdims=True)
  pred[qr]+=np.sum((corrected_train[tr]-pc_train[tr])[idx]*w[:,:,None],axis=1)
 return pred
def assign_labels(train,labels,query,k=15):
 from sklearn.neighbors import KNeighborsClassifier
 model=KNeighborsClassifier(n_neighbors=min(k,len(train)),weights='distance',n_jobs=4);model.fit(train,labels);return model.predict(query)
