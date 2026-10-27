from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,json,hashlib,datetime,gc,re
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase119_Bone_Lung_functional_subgroups_20261024'
for d in ['code','tables','objects','resources','logs','receipts','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(n,x):
 p=O/n;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8');t.replace(p)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def event(stage,status,**kw):
 d=dict(utc=now(),stage=stage,status=status,**kw)
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(d,ensure_ascii=False,default=str)+'\n')
 dump('STATUS.json',d);print(json.dumps(d,ensure_ascii=False,default=str),flush=True)
def load():return sparse.load_npz(O/'objects/log1p_cp10k.npz'),read('cohort_registry.tsv.gz'),read('genes.tsv').symbol.to_numpy()
def encode(site,i):return site+'_'+str(int(i)).zfill(2)
def technical(genes):return np.array([bool(re.match(r'^(MT-|RPL|RPS)',g)) or g in {'XIST','TSIX','MALAT1','DDX3Y','KDM5D','UTY','EIF1AY','RPS4Y1','RPS4Y2','ZFY','USP9Y'} for g in genes])
def fit_rep(X,m,genes,fitrows,transformrows,prefix,n_hvg=3000,npc=50):
 import scanpy as sc,anndata as ad,harmonypy as hm
 from sklearn.decomposition import PCA
 from sklearn.neighbors import NearestNeighbors
 cache=O/'objects'/(prefix+'_representation.npz');metadata=O/'resources'/(prefix+'_representation.json')
 if cache.exists() and metadata.exists():
  zz=np.load(cache);assert np.array_equal(zz['fit_rows'],fitrows) and np.array_equal(zz['transform_rows'],transformrows);assert set(zz['features'])<=set(genes)
  event(prefix,'reused_completed_representation',fit_cells=len(fitrows));return zz['embedding'],zz['raw_pca'],pd.Index(transformrows).get_indexer(fitrows)
 allowed=np.flatnonzero(~technical(genes));source_counts=m.iloc[fitrows].source.value_counts();eligible_sources=set(source_counts[source_counts>=2].index);hvgrows=fitrows[m.iloc[fitrows].source.isin(eligible_sources).to_numpy()]
 if len(hvgrows)<2:raise ValueError('Insufficient cells to estimate variable genes')
 dump('resources/'+prefix+'_HVG_source_eligibility.json',dict(source_counts=source_counts.to_dict(),HVG_eligible_sources=sorted(eligible_sources),HVG_only_excluded_rows=len(fitrows)-len(hvgrows),all_cells_retained_for_PCA_graph=True,reason='At least2 cells are mathematically required for within-source variability; single-cell sources cannot define HVGs'))
 a=ad.AnnData(X[hvgrows][:,allowed].copy(),obs=m.iloc[hvgrows][['source']].copy(),var=pd.DataFrame(index=genes[allowed]))
 sc.pp.highly_variable_genes(a,n_top_genes=min(n_hvg,len(allowed)),flavor='seurat',batch_key='source');cols=allowed[a.var.highly_variable.to_numpy()];del a
 raw=X[fitrows][:,cols].toarray().astype(np.float32);mu=raw.mean(0);sd=np.maximum(raw.std(0,ddof=1),1e-8);scaled=np.minimum((raw-mu)/sd,10);pca=PCA(n_components=min(npc,len(fitrows)-1,len(cols)-1),svd_solver='arpack',random_state=20261024);pc=pca.fit_transform(scaled).astype(np.float32)
 z=np.minimum((X[transformrows][:,cols].toarray()-mu)/sd,10);qpc=pca.transform(z).astype(np.float32);del raw,scaled,z;gc.collect()
 fm=m.iloc[fitrows];vc=fm.source.value_counts();use_harmony=len(vc)>1 and vc.min()>=10;reason='single_or_sparse_source';corrected=pc.copy();metric=[]
 if use_harmony:
  h=hm.run_harmony(pc,fm[['source']].reset_index(drop=True),['source'],theta=1,max_iter_harmony=20,random_state=20261024,ncores=4,verbose=False);hh=np.asarray(h.Z_corr);hh=hh.T if hh.shape==pc.T.shape else hh;assert hh.shape==pc.shape
  counts=fm.truth_lineage.value_counts();eligible=set(counts[(counts>=20)&(counts.index!='Unassessed')].index);good=np.flatnonzero(fm.truth_lineage.isin(eligible).to_numpy())
  if len(eligible)>=2:
   for name,emb in [('raw',pc),('harmony',hh)]:
    nm=min(15,len(good)-1);nn=NearestNeighbors(n_neighbors=nm+1,n_jobs=4).fit(emb[good]).kneighbors(emb[good],return_distance=False)[:,1:];lab=fm.iloc[good].truth_lineage.to_numpy();ss=fm.iloc[good].source.to_numpy();d=fm.iloc[good][['source','donor','truth_lineage']].copy();d['purity']=(lab[nn]==lab[:,None]).mean(1);d['cross_source_same_lineage']=((ss[nn]!=ss[:,None])&(lab[nn]==lab[:,None])).mean(1)
    k=d.groupby(['source','donor','truth_lineage'])[['purity','cross_source_same_lineage']].mean().groupby('truth_lineage').mean().reset_index();k['representation']=name;metric.append(k)
   md=pd.concat(metric);a=md[md.representation.eq('raw')].set_index('truth_lineage');b=md[md.representation.eq('harmony')].set_index('truth_lineage');loss=a.purity-b.purity;accepted=loss.mean()<=.05 and loss.max()<=.10;reason='passed_known_lineage_retention' if accepted else 'failed_known_lineage_retention_raw_fallback';save(md,prefix+'_representation_diagnostic.tsv')
  else:accepted=True;reason='insufficient_known_lineages_for_retention_test'
  if accepted:corrected=hh.astype(np.float32)
  else:use_harmony=False
 result=qpc.copy()
 if use_harmony:
  for src in m.iloc[transformrows].source.unique():
   tr=np.flatnonzero(fm.source.eq(src).to_numpy());qr=np.flatnonzero(m.iloc[transformrows].source.eq(src).to_numpy())
   if not len(tr):continue
   dist,ix=NearestNeighbors(n_neighbors=min(30,len(tr)),n_jobs=4).fit(pc[tr]).kneighbors(qpc[qr]);w=1/np.maximum(dist,1e-4);w/=w.sum(1,keepdims=True);result[qr]+=np.sum((corrected[tr]-pc[tr])[ix]*w[:,:,None],axis=1)
 # Exact training embeddings override their interpolated values.
 pos=pd.Index(transformrows).get_indexer(fitrows);assert (pos>=0).all();result[pos]=corrected
 save(pd.DataFrame({'symbol':genes[cols]}),prefix+'_HVG.tsv');dump('resources/'+prefix+'_representation.json',dict(fit_cells=len(fitrows),transform_cells=len(transformrows),genes=len(cols),dimensions=pc.shape[1],source_counts=vc.to_dict(),harmony=bool(use_harmony),diagnostic=reason,query_projection='Same-source nearest-neighbor learned displacement; approximation, not native Harmony query API'))
 np.savez_compressed(O/'objects'/(prefix+'_representation.npz'),fit_rows=fitrows,transform_rows=transformrows,embedding=result,raw_pca=qpc,features=genes[cols],mean=mu,sd=sd,components=pca.components_,pca_mean=pca.mean_)
 return result,qpc,pos
def graph_and_labels(z,res,seed=20261024,graph=None):
 import scanpy as sc,anndata as ad
 a=ad.AnnData(np.zeros((len(z),0),dtype=np.float32))
 if graph is None:
  a.obsm['rep']=z;sc.pp.neighbors(a,n_neighbors=min(15,len(z)-1),use_rep='rep',random_state=20261024);graph=a.obsp['connectivities'].tocsr()
 sc.tl.leiden(a,adjacency=graph,resolution=res,random_state=seed,flavor='leidenalg',key_added='cluster',n_iterations=-1);return a.obs.cluster.astype(str).to_numpy(),graph
def predict(z,labels,q):
 from sklearn.neighbors import KNeighborsClassifier
 k=KNeighborsClassifier(n_neighbors=min(15,len(z)),weights='distance',n_jobs=4).fit(z,labels);return k.predict(q)
def best_jaccard(mask,other):return float(max(np.sum(mask&(other==l))/np.sum(mask|(other==l)) for l in set(other)))
