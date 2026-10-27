from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,shutil
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score,silhouette_score
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase9_joint_20260919'
for d in ['tables','figures','objects','resources','delivery']:(O/d).mkdir(exist_ok=True)
SEED=20260919
protocol={'scope':'joint analysis of all previously accepted strict CAF candidates, not a forced reference-paper four-class annotation','expected_cells':4300,'primary_batch':'study_technology','primary_Harmony_theta':2,'theta_sensitivity':1,'resolution_primary':.5,'resolution_sensitivity':[.3,.8],'HVG':2000,'HVG_batch_key':'study_technology','PCs':30,'neighbors':15,'seed':SEED,'feature_alignment':'stable Ensembl IDs; published frozen gene mapping for symbol-only GSE250015; integration restricted to measured common features; no zero imputation for missing features','normalization':'log1p(CP10k) using full original library counts BEFORE common-gene selection','sensitivities':['uncorrected PCA','study_technology theta1','dataset-only theta2','sample-only theta2','no FLEX','exclude unconfirmed donor and GSE250015','DecontX available subset'],'biological_inference':'raw count pseudobulk, within each dataset; unknown donor excluded and same-patient samples aggregated; no pooled cross-cohort independent N','identity_labels_unchanged':True}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
(O/'STATUS.json').write_text(json.dumps({'status':'JOINT_INTEGRATION_RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
main=P/'PROJECT_STATE.json';shutil.copy2(main,P/'PROJECT_STATE.before_joint_20260919.json')
state=json.loads(main.read_text());state.update(current_stage=O.name,dynamic_status_file=O.name+'/STATUS.json',status='JOINT_INTEGRATION_RUNNING',current_status='JOINT_INTEGRATION_RUNNING',next_action='Complete explicitly authorized joint batch integration and biological/technical sensitivity review')
main.write_text(json.dumps(state,indent=2,ensure_ascii=False),encoding='utf-8')
files=[P/'phase3_identity_20260914/identity/objects/caf_strict.h5ad',P/'phase7_review_states_20260917/identity/objects/GSE152048_caf_strict.h5ad']+sorted((P/'geo_identity_20260919/objects').glob('*.h5ad'))
source=[];manifest=[];featuremap=[];common=None;corcounts=[];corids=[]
for f in files:
 a=ad.read_h5ad(f)
 dataset=str(a.obs.dataset.iloc[0])
 if dataset=='GSE250015':
  gm=pd.read_csv(P/'geo_cnv_20260918/inputs/GSE250015__GSM7969982/genes.tsv',sep='\t').set_index('source_gene_id')
  ids=gm.ensembl_id.reindex(a.var_names).to_numpy();assert pd.notna(ids).all()
 else:ids=np.array([v.split('.')[0] for v in a.var_names])
 assert len(set(ids))==len(ids)
 measured=set(v for v in ids if v.startswith('ENSG'))
 common=measured if common is None else common&measured
 o=pd.DataFrame(index=a.obs_names)
 for col in ['sample','patient_id','site','dataset']:
  o[col]=a.obs[col].astype(str).to_numpy()
 o['sample_key']=o.dataset+'|'+o['sample']
 o['technology']=a.obs['technology_reported'].astype(str).to_numpy()
 if 'technology' in a.obs:o['technology']=a.obs.technology.astype(str).to_numpy()
 o['study_technology']=o.dataset+'|'+o.technology
 o['donor_known']=~o.patient_id.isin(['UNCONFIRMED','unknown',''])
 o['donor_record']=np.where(o.donor_known,o.dataset+'|'+o.patient_id,o.sample_key+'|UNRESOLVED')
 o['inferential_eligible']=o.donor_known&(o.dataset!='GSE250015')
 o['full_library_counts']=np.asarray(a.layers['counts'].sum(1)).ravel()
 o['n_genes']=a.obs.n_genes.to_numpy()
 o['source_file']=str(f)
 o['legacy_identity_rule']=dataset=='GSE270231'
 o['corrected_available']=False
 with open(f,'rb') as fh:dig=hashlib.file_digest(fh,'sha256').hexdigest()
 manifest.append({'path':str(f),'dataset':dataset,'cells':a.n_obs,'genes':a.n_vars,'sha256':dig})
 featuremap.extend({'source_file':str(f),'source_feature':v,'ensembl_id':e,'symbol':str(s)} for v,e,s in zip(a.var_names,ids,a.var.symbol))
 corrected=None
 if 'decontx_counts' in a.layers:
  corrected=a.layers['decontx_counts'].copy();o['corrected_available']=True
 elif dataset=='GSE152048':
  blocks=[];names=[]
  for sample in o['sample'].unique():
   b=ad.read_h5ad(P/'phase7_review_states_20260917/remediation_20260917/decontx_comparison'/(sample+'_strict_corrected.h5ad'))
   ix=pd.Index(b.var_names).get_indexer(a.var_names);assert (ix>=0).all()
   blocks.append(b.X[:,ix]);names.extend(b.obs_names)
  c=sp.vstack(blocks).tocsr();order=pd.Index(names).get_indexer(a.obs_names);assert (order>=0).all()
  corrected=c[order];o['corrected_available']=True
 source.append({'a':a,'ids':ids,'obs':o,'corrected':corrected})
pd.DataFrame(manifest).to_csv(O/'resources/source_manifest.tsv',sep='\t',index=False)
pd.DataFrame(featuremap).to_csv(O/'resources/feature_mapping.tsv.gz',sep='\t',index=False)
genes=sorted(common);assert len(genes)>10000
symbols={}
for z in source:
 for e,s in zip(z['ids'],z['a'].var.symbol):
  if e not in symbols or str(z['a'].obs.dataset.iloc[0])=='GSE152048':symbols[e]=str(s)
blocks=[];obs=[];logs=[];cor=[];ci=[];ct=[]
programs=json.loads((P/'phase8_original_plan_20260919/resources/programs_frozen.json').read_text())
programs={k:v for k,v in programs.items() if not k.startswith('Frozen_')}
programs['S100A4_marker']=['S100A4']
coverage=[]
for z in source:
 a=z['a'];o=z['obs'];ix=pd.Index(z['ids']).get_indexer(genes);assert (ix>=0).all()
 x=a.layers['counts'][:,ix].tocsr();blocks.append(x);obs.append(o)
 y=x.astype(float).multiply((1e4/np.maximum(o.full_library_counts.to_numpy(),1))[:,None]).tocsr();y.data=np.log1p(y.data);logs.append(y)
 # Program evidence uses each source's complete measured gene space, not only common HVGs.
 full=a.layers['counts'].tocsc();sy=a.var.symbol.astype(str).to_numpy()
 for name,gs in programs.items():
  vals=np.zeros(a.n_obs);det=np.zeros(a.n_obs,dtype=int);n=0
  for g in gs:
   gi=np.flatnonzero(sy==g)
   if not len(gi):continue
   n+=1;v=np.asarray(full[:,gi].sum(1)).ravel()
   vals+=np.log1p(v/np.maximum(o.full_library_counts.to_numpy(),1)*1e4);det+=(v>0)
  o['program_'+name]=vals/max(n,1);o['det_'+name]=det
  coverage.append({'dataset':o.dataset.iloc[0],'source_file':o.source_file.iloc[0],'program':name,'covered':n,'requested':len(gs)})
 if z['corrected'] is not None:
  cor.append(z['corrected'][:,ix]);ci.extend(o.index);ct.extend(np.asarray(z['corrected'].sum(1)).ravel())
raw=sp.vstack(blocks).tocsr();allobs=pd.concat(obs);xlog=sp.vstack(logs).tocsr()
assert len(allobs)==4300 and allobs.index.is_unique
a=ad.AnnData(xlog,obs=allobs.copy(),var=pd.DataFrame({'symbol':[symbols[g] for g in genes]},index=genes))
a.layers['counts']=raw
a.uns['counts_scope']='unaltered raw counts for shared measured Ensembl genes; full libraries and source paths retained in obs; original full objects remain unchanged'
pd.DataFrame(coverage).to_csv(O/'tables/program_coverage.tsv',sep='\t',index=False)
a.obs.groupby(['dataset','sample','patient_id','donor_known','technology','site'],observed=True).size().rename('cells').reset_index().to_csv(O/'tables/sample_manifest.tsv',sep='\t',index=False)
a.var.to_csv(O/'tables/common_genes.tsv',sep='\t')
metrics=[];support=[];program_summary=[]
def pca(b):
 work=b.copy();sc.pp.filter_genes(work,min_cells=3)
 eligible=~work.var.symbol.str.startswith(('MT-','RPL','RPS')).to_numpy()
 work=work[:,eligible].copy()
 sc.pp.highly_variable_genes(work,n_top_genes=min(2000,work.n_vars),flavor='seurat',batch_key='study_technology')
 hv=work[:,work.var.highly_variable].copy();sc.pp.scale(hv,max_value=10)
 sc.tl.pca(hv,n_comps=min(30,hv.n_obs-1,hv.n_vars-1),svd_solver='arpack',random_state=SEED)
 return hv.obsm['X_pca'].copy(),hv.var_names.tolist()
def harmony(base,b,key,theta):
 h=hm.run_harmony(base,b.obs,key,theta=theta,random_state=SEED,ncores=4,max_iter_harmony=30,verbose=True)
 v=np.asarray(h.Z_corr)
 if v.shape!=base.shape:v=v.T
 assert v.shape==base.shape and np.isfinite(v).all()
 return v
def cluster(b,rep,name,res=.5,reference=None,umap=False):
 b.obsm[name]=rep;sc.pp.neighbors(b,use_rep=name,n_neighbors=15,random_state=SEED)
 sc.tl.leiden(b,resolution=res,key_added=name,random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
 labs=b.obs[name].astype(str).to_numpy();batch=b.obs.study_technology.astype(str).to_numpy()
 nn=NearestNeighbors(n_neighbors=16).fit(rep).kneighbors(rep,return_distance=False)[:,1:]
 row={'analysis':name,'cells':b.n_obs,'clusters':len(set(labs)),'batch_silhouette':float(silhouette_score(rep,batch,sample_size=min(2000,b.n_obs),random_state=SEED)),'same_batch_neighbor_fraction':float((batch[nn]==batch[:,None]).mean()),'same_donor_record_neighbor_fraction':float((b.obs.donor_record.to_numpy()[nn]==b.obs.donor_record.to_numpy()[:,None]).mean())}
 if reference is not None:row['ARI_vs_main']=adjusted_rand_score(reference,labs)
 metrics.append(row)
 b.obs[['dataset','sample','patient_id','donor_record','donor_known','site',name]].to_csv(O/'tables'/(name+'_cell_labels.tsv.gz'),sep='\t',index_label='cell_id')
 for cl in sorted(set(labs),key=int):
  mask=labs==cl;df=b.obs.iloc[np.flatnonzero(mask)];known=df[df.donor_known]
  counts=known.groupby(['dataset','donor_record'],observed=True).size()
  record={'analysis':name,'cluster':cl,'cells':int(mask.sum()),'datasets':int(df.dataset.nunique()),'max_dataset_fraction':float(df.dataset.value_counts().max()/len(df)),'known_donor_records_ge10':int((counts>=10).sum()),'known_donor_records_ge20':int((counts>=20).sum()),'max_known_donor_fraction':float(counts.max()/len(df)) if len(counts) else np.nan,'unknown_donor_cells':int((~df.donor_known).sum()),'lung_cells':int(df.site.str.lower().eq('lung').sum())}
  support.append(record)
  for col in [c for c in df if c.startswith('program_')]:
   program_summary.append({'analysis':name,'cluster':cl,'program':col[8:],'mean_score':float(df[col].mean()),'median_score':float(df[col].median())})
 if umap:
  sc.tl.umap(b,random_state=SEED);b.obsm['X_umap_'+name]=b.obsm['X_umap'].copy()
 return labs
base,hvg=pca(a);a.obsm['X_pca_joint']=base
(O/'resources/HVG_primary.json').write_text(json.dumps(hvg))
cluster(a,base,'uncorrected',umap=True)
corrected=harmony(base,a,'study_technology',2)
main_labels=cluster(a,corrected,'joint_main',umap=True)
main_graph=a.obsp['connectivities'].copy();main_dist=a.obsp['distances'].copy();main_neighbors=a.uns['neighbors'].copy()
for resolution in [.3,.8]:cluster(a,corrected,'joint_resolution_'+str(resolution),res=resolution,reference=main_labels)
for key,theta,name in [('study_technology',1,'joint_theta1'),('dataset',2,'joint_dataset_only'),('sample_key',2,'joint_sample_batch')]:
 v=harmony(base,a,key,theta);cluster(a,v,name,reference=main_labels)
a.obs['joint_cluster']=a.obs.joint_main.astype(str)
for label,mask in [('no_FLEX',~a.obs.technology.str.contains('FLEX')),('mapped_without_GSE250015',a.obs.donor_known&(a.obs.dataset!='GSE250015'))]:
 b=a[mask].copy();bp,bhv=pca(b);rep=harmony(bp,b,'study_technology',2)
 cluster(b,rep,label,reference=a.obs.loc[b.obs_names,'joint_main'].astype(str).to_numpy())
 b.write_h5ad(O/'objects'/(label+'.h5ad'),compression='gzip')
# Corrected-count sensitivity excludes the original 188 cells lacking matched DecontX.
cm=sp.vstack(cor).tocsr();ct=np.asarray(ct);order=pd.Index(ci).get_indexer(a.obs_names[a.obs.corrected_available]);assert (order>=0).all()
b=a[a.obs.corrected_available].copy();cr=cm[order]
b.layers['decontx_counts']=cr;b.X=cr.astype(float).multiply((1e4/np.maximum(ct[order],1))[:,None]).tocsr();b.X.data=np.log1p(b.X.data)
bp,bhv=pca(b);rep=harmony(bp,b,'study_technology',2)
cluster(b,rep,'decontx_available',reference=a.obs.loc[b.obs_names,'joint_main'].astype(str).to_numpy())
b.write_h5ad(O/'objects/decontx_available.h5ad',compression='gzip')
a.obsp['connectivities']=main_graph;a.obsp['distances']=main_dist;a.uns['neighbors']=main_neighbors;a.obsm['X_umap']=a.obsm['X_umap_joint_main'].copy()
a.write_h5ad(O/'objects/joint_CAF_integrated.h5ad',compression='gzip')
check=ad.read_h5ad(O/'objects/joint_CAF_integrated.h5ad')
assert check.obs_names.equals(a.obs_names) and (check.layers['counts']!=raw).nnz==0
pd.DataFrame(metrics).to_csv(O/'tables/integration_metrics.tsv',sep='\t',index=False)
pd.DataFrame(support).to_csv(O/'tables/cluster_support.tsv',sep='\t',index=False)
pd.DataFrame(program_summary).to_csv(O/'tables/cluster_programs.tsv',sep='\t',index=False)
(O/'resources/programs.json').write_text(json.dumps(programs,indent=2))
summary={'status':'JOINT_INTEGRATION_COMPUTED_REVIEW_PENDING','cells':a.n_obs,'shared_genes':a.n_vars,'datasets':a.obs.dataset.nunique(),'study_technology_batches':a.obs.study_technology.nunique(),'sample_records':a.obs.sample_key.nunique(),'known_donor_records_not_pooled_independent_N':a.obs.loc[a.obs.donor_known,'donor_record'].nunique(),'unknown_donor_cells':int((~a.obs.donor_known).sum()),'main_clusters':a.obs.joint_main.nunique(),'original_counts_unchanged':True}
(O/'STATUS.json').write_text(json.dumps(summary,indent=2))
(O/'INTEGRATION_COMPLETE').write_text('JOINT_INTEGRATION_COMPUTED_REVIEW_PENDING\n')
print(json.dumps(summary),flush=True);print(pd.DataFrame(metrics).to_string(index=False),flush=True)
