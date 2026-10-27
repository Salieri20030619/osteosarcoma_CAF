from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,fcntl
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase27_identity_review_20260926';O=P/'phase28b_CAF_refined_20260926'
for d in ['code','tables','objects','figures','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(p,v):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str));t.replace(p)
def status(s,**kw):writej(O/'STATUS.json',{'status':s,'at':now(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def clean(a):
 for frame in [a.obs,a.var]:
  for c in frame:
   if frame[c].dtype==object or isinstance(frame[c].dtype,pd.StringDtype):frame[c]=frame[c].fillna('').astype(str).astype(object)
  frame.index=pd.Index(frame.index.astype(str).to_numpy(dtype=object))
def model(a):
 status('NORMALIZATION_PCA',cells=a.n_obs,genes=a.n_vars)
 a.layers['counts']=a.X.copy()
 a.X=a.X.astype(np.float32).multiply((1e4/a.obs.full_library_counts.to_numpy())[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
 sc.pp.highly_variable_genes(a,n_top_genes=min(2000,a.n_vars),flavor='seurat',batch_key='study_technology')
 w=a[:,a.var.highly_variable].copy();sc.pp.scale(w,max_value=10)
 sc.tl.pca(w,n_comps=min(50,w.n_obs-1,w.n_vars-1),svd_solver='arpack',random_state=20260926)
 z=w.obsm['X_pca'];a.obsm['X_pca']=z
 status('HARMONY_NEIGHBORS',cells=a.n_obs)
 if a.obs.study_technology.nunique()>1:
  h=hm.run_harmony(z,a.obs,'study_technology',random_state=20260926,ncores=4,verbose=False)
  h=np.asarray(h.Z_corr);h=h.T if h.shape!=z.shape else h;assert h.shape==z.shape and np.isfinite(h).all()
 else:h=z.copy()
 a.obsm['X_harmony']=h
 sc.pp.neighbors(a,use_rep='X_harmony',n_neighbors=min(15,a.n_obs-1),random_state=20260926)
 sc.tl.leiden(a,resolution=1,key_added='cluster',flavor='igraph',n_iterations=2,directed=False,random_state=20260926)
 sc.tl.leiden(a,resolution=.5,key_added='cluster_res05',flavor='igraph',n_iterations=2,directed=False,random_state=20260926)
 sc.tl.umap(a,random_state=20260926)
 return a
def profiles(a,prefix='CAF'):
 clusters=sorted(a.obs.cluster.unique(),key=lambda x:int(x));rows=[];markers=[];programs=[]
 names=[c for c in a.obs if c.startswith(('lineage_','program_'))]
 for cl in clusters:
  mask=(a.obs.cluster==cl).to_numpy();g=a.obs.loc[mask];dn=g.groupby('donor_record',observed=True).size()
  row=dict(cluster=cl,cells=int(mask.sum()),sample_records=g['sample'].nunique(),donor_records=len(dn),datasets=g.dataset.nunique(),largest_donor_fraction=dn.max()/mask.sum(),top_donor=dn.idxmax(),dual_supported=int(g.CNV_evidence.eq('concordant_no_alteration').sum()),single_supported=int(g.CNV_evidence.eq('single_method_no_alteration').sum()))
  # Donor-equal panel means describe each cluster; no cell-level significance tests.
  for n in names:row[n]=float(g.groupby('donor_record',observed=True)[n].mean().mean())
  row['provider_labels']=';'.join(f'{k}:{v}' for k,v in g.provider_annotation.value_counts().head(6).items())
  row['historical_labels']=';'.join(f'{k}:{v}' for k,v in g.historical_annotation.value_counts().head(6).items())
  rows.append(row)
  for keys,q in g.groupby(['dataset','donor_record','site','clinical_group'],observed=True):
   programs.append(dict(cluster=cl,dataset=keys[0],donor_record=keys[1],site=keys[2],clinical_group=keys[3],cells=len(q),**{n:float(q[n].mean()) for n in names}))
  mean=np.asarray(a.X[mask].mean(0)).ravel();other=np.asarray(a.X[~mask].mean(0)).ravel() if (~mask).any() else mean
  frac=a.layers['counts'][mask].getnnz(0)/mask.sum();fraco=a.layers['counts'][~mask].getnnz(0)/(~mask).sum() if (~mask).any() else frac
  ix=np.argsort(mean-other)[::-1][:100]
  for j in ix:markers.append(dict(cluster=cl,ensembl=a.var_names[j],symbol=a.var.symbol.iloc[j],mean_logCP10k=mean[j],rest_mean_logCP10k=other[j],delta_logCP10k=mean[j]-other[j],detection_fraction=frac[j],rest_detection_fraction=fraco[j],ranking='descriptive_cell_expression_not_inferential'))
 save(pd.DataFrame(rows),prefix+'_cluster_review.tsv');save(pd.DataFrame(markers),prefix+'_descriptive_markers.tsv');save(pd.DataFrame(programs),prefix+'_donor_programs.tsv')
 return pd.DataFrame(rows)
def plots(a,profile):
 fig,axes=plt.subplots(1,3,figsize=(17,5),layout='constrained')
 for ax,col in zip(axes,['cluster','dataset','CNV_evidence']):
  for label,g in a.obs.groupby(col,observed=True):
   ix=a.obs.index.get_indexer(g.index);ax.scatter(a.obsm['X_umap'][ix,0],a.obsm['X_umap'][ix,1],s=2,label=str(label),rasterized=True)
  ax.set_title(col);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(markerscale=3,fontsize=6,loc='best',frameon=False)
 fig.savefig(O/'figures/CAF_umap.png',dpi=170);plt.close(fig)
 cols=[c for c in profile if c.startswith(('lineage_','program_'))]
 vals=profile[cols].to_numpy();zs=(vals-vals.mean(0))/np.maximum(vals.std(0),1e-8)
 fig,ax=plt.subplots(figsize=(14,max(5,len(profile)*.35)),layout='constrained');im=ax.imshow(zs,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
 ax.set_yticks(range(len(profile)),['C'+str(c) for c in profile.cluster]);ax.set_xticks(range(len(cols)),[c.replace('lineage_','L:').replace('program_','P:') for c in cols],rotation=70,ha='right');fig.colorbar(im,ax=ax,label='Across-cluster standardized donor-mean score');ax.set_title('CAF candidate cluster review')
 fig.savefig(O/'figures/CAF_program_review.png',dpi=160);plt.close(fig)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);rc=1
try:
 assert (B/'REVIEW_ACCEPTED').exists() and (B/'acceptance.exit_code').read_text().strip()=='0'
 old=P/'phase28_CAF_clustering_20260926'
 a=ad.read_h5ad(old/'objects/CAF_integrated.h5ad')
 annotation=pd.read_csv(old/'tables/CAF_cluster_annotation_reviewed.tsv',sep='\t')
 retain=annotation.loc[annotation.identity_review.eq('retain_CAF_candidate'),'cluster'].astype(str)
 a=a[a.obs.cluster.astype(str).isin(retain)].copy()
 a.obs['identity_review_parent_cluster']=a.obs.cluster.astype(str)
 a.obs['identity_overlap_flag']=a.obs.cluster.astype(str).eq('14')
 a.X=a.layers['counts'].copy()
 cm=pd.read_csv(B/'tables/clinical_mapping_v3.tsv',sep='\t').fillna('').set_index('sample')
 for c in ['clinical_group','group_evidence']:a.obs[c]=a.obs['sample'].map(cm[c])
 assert a.n_obs==1454
 assert a.obs.CNV_evidence.isin(['single_method_no_alteration','concordant_no_alteration']).all()
 writej(O/'PROTOCOL.json',{'input':'1454 fibroblast candidates retained after 28-cluster biological identity review; parent C14 osteogenic overlap retained with sensitivity flag','cell_entry_filters_added':False,'normalization':'log1p CP10k using each full original library denominator','HVG':2000,'PCA':50,'batch':'study_technology','Harmony':'installed defaults, seed20260926, 4cores','Leiden':'resolution1 primary, resolution0.5 boundary sensitivity, igraph 2iterations undirected','biological_P_values':'none in marker ranking; descriptive cluster expression only','next':'manual cluster-level CAF typing from customer panels and lineage review, then myCAF internal refit','original_data_preserved':True})
 a=model(a);profile=profiles(a);plots(a,profile)
 clean(a);a.write_h5ad(O/'objects/CAF_integrated.h5ad',compression='gzip');a.obs.to_csv(O/'tables/CAF_cell_labels.tsv.gz',sep='\t',index_label='cell_id')
 save(a.obs.groupby(['sample','dataset','donor_record','clinical_group','cluster'],observed=True).size().reset_index(name='cells'),'CAF_sample_cluster_counts.tsv')
 qa={'at':now(),'cells':a.n_obs,'genes':a.n_vars,'clusters':a.obs.cluster.nunique(),'clusters_resolution05':a.obs.cluster_res05.nunique(),'sample_records':a.obs['sample'].nunique(),'donor_records_not_crossstudy_independent_N':a.obs.donor_record.nunique(),'finite_embedding':bool(np.isfinite(a.obsm['X_umap']).all()),'CNV_conflict_or_unassessed_in_main':0,'myCAF_reclustering':'pending cluster annotation review'}
 writej(O/'delivery/CLUSTER_QA.json',qa);(O/'CLUSTERING_COMPLETE').write_text(now()+'\n');status('CAF_CLUSTERING_COMPLETE_ANNOTATION_REVIEW_REQUIRED',**qa);rc=0
except Exception as e:status('FAILED',error=str(e));traceback.print_exc()
finally:(O/'exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)
