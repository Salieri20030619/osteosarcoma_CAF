from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,warnings
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse as sp
from sklearn.metrics import adjusted_rand_score,silhouette_score
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(4);warnings.simplefilter('ignore',pd.errors.PerformanceWarning)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';assert (O/'QUALITY_GUARD_COMPLETE').exists();SEED=20260919
protocol={'primary':'unified_strict, raw-count logCP10k, HVG2000 PC30, Harmony study_technology theta2 dynamic lambda, k15 Leiden0.5','interpretation':'predeclared working baseline, not selected by maximum mixing or significance','HVG_PC_grid':{'HVG':[1500,2000,3000],'PC':[20,30,40]},'sensitivity':['uncorrected','theta1','sample-level batch','fixed_lambda1','DecontX counts','exclude FLEX','balanced PCA fit cap200 per donor record','osteo threshold .5','all competition thresholds .5','leave-one-known-donor-record-out refit'],'batch_metrics':'silhouette, same-batch neighbors with composition expectation; not sole tuning targets','biological_check':'raw program neighbor distance and subset label stability are descriptive, not external validation','donor_policy':'Unknown donor retained in embedding but not grouped as biological replicates; crosscohort records not pooled N'}
(O/'resources/batch_protocol.json').write_text(json.dumps(protocol,indent=2))
man=pd.read_csv(O/'resources/object_manifest.tsv',sep='\t');sources=[];common=None;symbols={}
for row in man.to_dict('records'):
 b=ad.read_h5ad(row['path']);assert not b.obs.independent_doublet_excluded.any();mask=b.obs[['unified_strict','osteo_sensitive','competitor_sensitive']].any(axis=1).to_numpy()
 if not mask.any():continue
 b=b[mask].copy();ids=set(g for g in b.var_names if g.startswith('ENSG'));common=ids if common is None else common&ids
 for g,s in zip(b.var_names,b.var.symbol.astype(str)):
  if g not in symbols or row['dataset']=='GSE152048':symbols[g]=s
 sources.append(b)
genes=sorted(common);assert len(genes)>10000
raw=[];corrected=[];obs=[]
for b in sources:
 ix=b.var_names.get_indexer(genes);assert (ix>=0).all();raw.append(b.layers['counts'][:,ix]);corrected.append(b.layers['decontx_counts'][:,ix]);obs.append(b.obs.copy())
x=sp.vstack(raw).tocsr();cx=sp.vstack(corrected).tocsr();ob=pd.concat(obs);assert ob.index.is_unique
var=pd.DataFrame({'symbol':[symbols[g] for g in genes]},index=genes)
def normalize(x,den):
 y=x.astype(float).multiply((1e4/np.maximum(np.asarray(den),1))[:,None]).tocsr();y.data=np.log1p(y.data);return y
allobj=ad.AnnData(normalize(x,ob.full_library_counts),obs=ob,var=var);allobj.layers['counts']=x;allobj.layers['decontx_counts']=cx;allobj.uns['count_scope']='unaltered raw/corrected source counts for common measured Ensembl genes; full libraries in obs'
allobj.write_h5ad(O/'objects/unified_candidate_union.h5ad',compression='gzip')
a=allobj[allobj.obs.unified_strict].copy();assert a.n_obs>100
metrics=[];support=[];overlap=[];leaveout=[]

# PCA and scaling refit per subset; balanced option learns HVGs/scaling/PCA on capped donor records then projects all cells.
def pca(b,hvg=2000,npc=30,balanced=False):
 w=b[:,~b.var.symbol.str.startswith(('MT-','RPL','RPS')).to_numpy()].copy();sc.pp.filter_genes(w,min_cells=3)
 if not balanced:
  sc.pp.highly_variable_genes(w,n_top_genes=min(hvg,w.n_vars),flavor='seurat',batch_key='study_technology');v=w[:,w.var.highly_variable].copy();sc.pp.scale(v,max_value=10);sc.tl.pca(v,n_comps=min(npc,v.n_obs-1,v.n_vars-1),svd_solver='arpack',random_state=SEED);return v.obsm['X_pca'].copy(),v.var_names.tolist()
 from sklearn.preprocessing import StandardScaler
 from sklearn.decomposition import PCA
 rng=np.random.default_rng(SEED);indices=[]
 for _,ii in w.obs.groupby('donor_record',observed=True).indices.items():indices.extend(rng.choice(ii,min(200,len(ii)),replace=False))
 indices=np.array(sorted(indices));t=w[indices].copy();sc.pp.highly_variable_genes(t,n_top_genes=min(hvg,t.n_vars),flavor='seurat',batch_key='study_technology');gv=t.var_names[t.var.highly_variable];v=w[:,gv].X.toarray();scaler=StandardScaler().fit(v[indices]);z=np.clip(scaler.transform(v),-10,10);model=PCA(n_components=min(npc,len(indices)-1,len(gv)-1),svd_solver='arpack',random_state=SEED).fit(z[indices]);return model.transform(z),gv.tolist()

def harmonize(z,b,key='study_technology',theta=2,lamb=None):
 h=hm.run_harmony(z,b.obs,key,theta=theta,lamb=lamb,random_state=SEED,ncores=4,max_iter_harmony=30,verbose=False);v=np.asarray(h.Z_corr);v=v.T if v.shape!=z.shape else v;assert v.shape==z.shape and np.isfinite(v).all();return v

def evaluate(b,z,name,reference=None,umap=False,seed=SEED):
 b.obsm[name]=z;sc.pp.neighbors(b,use_rep=name,n_neighbors=15,random_state=SEED);sc.tl.leiden(b,resolution=.5,key_added=name,random_state=seed,flavor='igraph',n_iterations=2,directed=False);labels=b.obs[name].astype(str);batch=b.obs.study_technology.astype(str).to_numpy();nn=NearestNeighbors(n_neighbors=16).fit(z).kneighbors(z,return_distance=False)[:,1:]
 probs=pd.Series(batch).value_counts(normalize=True).to_numpy();expected=float((probs**2).sum());same=float((batch[nn]==batch[:,None]).mean())
 programs=[g for g in ['Fibroblast_core','Contractile','Matrix_remodeling','Inflammatory','Antigen_presentation','Osteogenic_full_customer','Calcium_regulation_full_customer','Cell_cycle'] if 'program_'+g in b.obs];vals=b.obs[['program_'+g for g in programs]].to_numpy(float);std=vals.std(0);vals=(vals-vals.mean(0))/np.maximum(std,1e-9);dist=float(np.abs(vals[nn]-vals[:,None,:]).mean())
 rec={'analysis':name,'cells':b.n_obs,'clusters':labels.nunique(),'batch_silhouette':float(silhouette_score(z,batch,sample_size=min(1500,b.n_obs),random_state=SEED)) if len(set(batch))>1 else np.nan,'same_batch_neighbors':same,'composition_expected_same_batch':expected,'same_batch_excess':same-expected,'raw_program_neighbor_distance':dist,'ARI_vs_unified_primary':adjusted_rand_score(reference,labels) if reference is not None else np.nan};metrics.append(rec)
 b.obs[['dataset','sample','patient_id','donor_record','donor_known','technology','site',name]].to_csv(O/'tables'/(name+'_labels.tsv.gz'),sep='\t',index_label='cell_id')
 for cl,g in b.obs.groupby(name,observed=True):
  counts=g[g.donor_known].groupby('donor_record',observed=True).size();support.append({'analysis':name,'cluster':str(cl),'cells':len(g),'datasets':g.dataset.nunique(),'donor_records_ge20':int((counts>=20).sum()),'max_known_donor_fraction':float(counts.max()/len(g)) if len(counts) else np.nan,'unknown_donor_cells':int((~g.donor_known).sum())})
 if reference is not None:
  rr=np.asarray(reference).astype(str);ll=labels.to_numpy()
  for cl in sorted(set(rr),key=int):
   u=rr==cl;best=max((np.sum(u&(ll==v))/np.sum(u|(ll==v)),v) for v in set(ll));overlap.append({'analysis':name,'primary_cluster':cl,'primary_cells_shared':int(u.sum()),'best_cluster':best[1],'Jaccard':float(best[0])})
 if umap:sc.tl.umap(b,random_state=SEED);b.obsm['X_umap_'+name]=b.obsm['X_umap'].copy()
 print('BATCH_RUN',json.dumps(rec),flush=True);return labels.to_numpy()

z,hv=pca(a);a.obsm['X_pca']=z;(O/'resources/main_HVG.json').write_text(json.dumps(hv))
evaluate(a,z,'unified_uncorrected',umap=True);v=harmonize(z,a);main=evaluate(a,v,'unified_primary',umap=True);a.obs['unified_cluster']=main
main_graph=a.obsp['connectivities'].copy();main_dist=a.obsp['distances'].copy();main_neighbors=a.uns['neighbors'].copy()
for h in [1500,2000,3000]:
 for n in [20,30,40]:
  if h==2000 and n==30:continue
  zz,_=pca(a,h,n);evaluate(a,harmonize(zz,a),f'hvg{h}_pc{n}',reference=main)
for key,theta,lam,name in [('study_technology',1,None,'theta1'),('sample_key',2,None,'sample_batch'),('study_technology',2,1,'fixed_lambda1')]:evaluate(a,harmonize(z,a,key=key,theta=theta,lamb=lam),name,reference=main)
zz,_=pca(a,balanced=True);evaluate(a,harmonize(zz,a),'balanced_PCA_fit',reference=main)
evaluate(a,v,'seed_sensitivity',reference=main,seed=20260920)
for name,mask,correction in [('no_FLEX',~a.obs.technology.str.contains('FLEX'),False),('corrected_counts',np.ones(a.n_obs,bool),True)]:
 b=a[mask].copy()
 if correction:b.X=normalize(b.layers['decontx_counts'],b.obs.corrected_library_counts)
 zz,_=pca(b);evaluate(b,harmonize(zz,b),name,reference=a.obs.loc[b.obs_names,'unified_cluster'].to_numpy(),umap=True);b.write_h5ad(O/'objects'/(name+'_integrated.h5ad'),compression='gzip')
for name in ['osteo_sensitive','competitor_sensitive']:
 b=allobj[allobj.obs[name]].copy();zz,_=pca(b);labs=evaluate(b,harmonize(zz,b),name,umap=True)
 shared=a.obs_names.intersection(b.obs_names);rr=a.obs.loc[shared,'unified_cluster'].astype(str);ll=b.obs.loc[shared,name].astype(str)
 for cl in sorted(rr.unique(),key=int):
  u=rr==cl;best=max((int((u&(ll==v)).sum())/int((u|(ll==v)).sum()),v) for v in ll.unique());overlap.append({'analysis':name,'primary_cluster':cl,'primary_cells_shared':int(u.sum()),'best_cluster':best[1],'Jaccard':float(best[0])})
 metrics[-1]['ARI_vs_unified_primary']=adjusted_rand_score(rr,ll);b.write_h5ad(O/'objects'/(name+'_integrated.h5ad'),compression='gzip')
# Leave known donor records out; this tests dependence on a record, not clinical independence or predictive validation.
for donor in sorted(a.obs.loc[a.obs.donor_known,'donor_record'].unique()):
 b=a[a.obs.donor_record!=donor].copy()
 if b.n_obs<100:continue
 zz,_=pca(b);name='leaveout_'+str(len(leaveout));ref=a.obs.loc[b.obs_names,'unified_cluster'].to_numpy();before=len(overlap);evaluate(b,harmonize(zz,b),name,reference=ref)
 for row in overlap[before:]:row['omitted_donor_record']=donor
 leaveout.append({'analysis':name,'omitted_donor_record':donor,'omitted_cells':int((a.obs.donor_record==donor).sum()),'remaining_cells':b.n_obs,'ARI':metrics[-1]['ARI_vs_unified_primary']})
a.obsp['connectivities']=main_graph;a.obsp['distances']=main_dist;a.uns['neighbors']=main_neighbors;a.obsm['X_umap']=a.obsm['X_umap_unified_primary'].copy();a.write_h5ad(O/'objects/unified_strict_integrated.h5ad',compression='gzip')
pd.DataFrame(metrics).to_csv(O/'tables/batch_model_comparison.tsv',sep='\t',index=False);pd.DataFrame(support).to_csv(O/'tables/cluster_support.tsv',sep='\t',index=False);pd.DataFrame(overlap).to_csv(O/'tables/model_cluster_overlap.tsv',sep='\t',index=False);pd.DataFrame(leaveout).to_csv(O/'tables/leaveout_summary.tsv',sep='\t',index=False)
check=ad.read_h5ad(O/'objects/unified_strict_integrated.h5ad');assert check.obs_names.equals(a.obs_names) and (check.layers['counts']!=a.layers['counts']).nnz==0
status=json.loads((O/'STATUS.json').read_text());status.update(status='UNIFIED_GATES_AND_BATCH_MODELS_COMPLETE_REVIEW_PENDING',shared_genes=len(genes),main_clusters=len(set(main)),batch_models=len(metrics),main_datasets=int(a.obs.dataset.nunique()),main_samples=int(a.obs.sample_key.nunique()),main_cells=a.n_obs);(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'BATCH_MODELS_COMPLETE').write_text('BATCH_MODELS_COMPLETE\n');print(json.dumps(status),flush=True)