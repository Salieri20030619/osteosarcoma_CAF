from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,gzip,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy.io import mmread
from scipy import sparse
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase32_lung_sensitivity_20260926';M=P/'phase30_two_layer_classification_20260926';Q=P/'phase31_customer_completion_20260926';B=P/'phase26_client_reanalysis_20260922'
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):
 t=f.with_name(f.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'SENSITIVITY_STATUS.json',dict(status=s,at=now(),**kw));print(s,kw,flush=True)
def clean(a):
 for fr in [a.obs,a.var]:
  for c in fr:
   if fr[c].dtype==object or isinstance(fr[c].dtype,pd.StringDtype):fr[c]=fr[c].fillna('').astype(str).astype(object)
  fr.index=pd.Index(fr.index.astype(str).to_numpy(dtype=object))
scenarios={'decontx':[], 'raw_adjust_cycle':['S_score','G2M_score'],'raw_adjust_stress':['stress512_signed'],'raw_adjust_stress_cycle':['S_score','G2M_score','stress512_signed'],'decontx_adjust_stress_cycle':['S_score','G2M_score','stress512_signed']}
protocol={'at':now(),'cells_fixed':{'supported':1224,'expanded':8174,'supported_myCAF':90,'expanded_myCAF':225},'scenarios':scenarios,'baseline':'existing phase30 per-layer refit','model_parameters':'same HVG2000 seurat by study_technology, PCA<=50, Harmony defaults, neighbors<=15, Leiden1 and .5, seed20260926','covariates':'frozen full512 and formal S/G2M scores from original full source count space; continuous adjustment, no threshold-based cell removal','raw_counts_retained':True,'corrected_counts_are_estimates_not_integer_UMIs':True,'regression':'Scanpy regress_out continuous covariates on selected variable log-normalized genes, before scaling; raw marker expression is not replaced by residuals','acceptance':'same cell IDs, finite embeddings, both resolutions, original-to-refit contingency, per-class dispersion and ARI/NMI; no invented pass fraction','not_validated_CAF_lineages':True}
writej(O/'resources/SENSITIVITY_PROTOCOL.json',protocol)
try:
 orig=ad.read_h5ad(M/'objects/expanded_integrated.h5ad');ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id').loc[orig.obs_names]
 assert orig.n_obs==8174
 cf=O/'objects/expanded_decontx_counts.h5ad'
 if not cf.exists():
  mats=[];ids=[];umis=[];coverage=[]
  for sm,sub in ev.groupby('sample',sort=False):
   status('ASSEMBLING_CORRECTED_COUNTS',sample=sm,assembled_cells=len(ids))
   genes=read(B/'inputs'/sm/'genes.tsv');cells=read(B/'inputs'/sm/'cells.tsv');cont=read(Q/'decontx'/sm/'cell_contamination.tsv').set_index('cell_id')
   assert cells.cell_id.tolist()==cont.index.tolist()
   gi=pd.Index(genes.ensembl).get_indexer(orig.var_names);ci=pd.Index(cells.cell_id).get_indexer(sub.index)
   assert (gi>=0).all() and (ci>=0).all()
   with gzip.open(Q/'decontx'/sm/'corrected_counts.mtx.gz','rb') as f:
    with threadpool_limits(limits=4):y=mmread(f).tocsr()
   assert y.shape==(len(genes),len(cells))
   mats.append(y[gi,:][:,ci].T.tocsr().astype(np.float32));ids.extend(sub.index);umis.extend(cont.loc[sub.index,'corrected_UMI'].tolist());coverage.append(dict(sample=sm,candidate_cells=len(sub),common_genes=len(gi)))
   del y;gc.collect()
  c=ad.AnnData(sparse.vstack(mats,format='csr'),obs=ev.loc[ids].copy(),var=orig.var.copy());c.obs['corrected_full_library_counts']=umis;c=c[orig.obs_names].copy()
  assert np.isfinite(c.X.data).all() and (c.X.data>=0).all()
  clean(c);c.write_h5ad(cf,compression='gzip');pd.DataFrame(coverage).to_csv(O/'tables/corrected_assembly_coverage.tsv',sep='\t',index=False)
  del mats;gc.collect()
 else:c=ad.read_h5ad(cf)
 scores=read(Q/'stress_cycle/full512_and_cell_cycle_scores.tsv.gz').set_index('cell_id')
 summaries=[];profiles=[];correspond=[];outcomes=[]
 for layer in ['supported','expanded','supported_myCAF','expanded_myCAF']:
  base=ad.read_h5ad(M/f'objects/{layer}_integrated.h5ad')
  labfile=f'{layer}_cell_labels.tsv.gz' if layer.endswith('myCAF') else f'{layer}_classified_cells.tsv.gz'
  labels=read(M/'tables'/labfile).set_index('cell_id').loc[base.obs_names]
  for scenario,covs in scenarios.items():
   tag=layer+'__'+scenario;done=O/'tables'/f'{tag}_QA.json'
   if done.exists():summaries.append(json.loads(done.read_text()));outcomes.append(dict(layer=layer,scenario=scenario,status='complete_cached'));continue
   try:
    status('FITTING',layer=layer,scenario=scenario,cells=base.n_obs,completed=len(outcomes))
    a=ad.AnnData((c[base.obs_names,base.var_names].X if scenario.startswith('decontx') else base.layers['counts']).copy(),obs=labels.copy(),var=base.var.copy())
    denom=c.obs.loc[a.obs_names,'corrected_full_library_counts'].to_numpy(float) if scenario.startswith('decontx') else a.obs.full_library_counts.to_numpy(float)
    a.obs['normalization_library_counts']=denom;a.obs['original_cluster']=base.obs.cluster.astype(str).to_numpy()
    a.obs['frozen_class_id']=labels.class_id.to_numpy()
    for col in ['S_score','G2M_score','stress512_signed']:a.obs[col]=scores.loc[a.obs_names,col].to_numpy(float)
    assert np.isfinite(a.obs[['S_score','G2M_score','stress512_signed']].to_numpy()).all()
    a.X=a.X.astype(np.float32).multiply((1e4/np.maximum(denom,1))[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
    sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='study_technology')
    w=a[:,a.var.highly_variable].copy()
    if covs:sc.pp.regress_out(w,keys=covs,n_jobs=4)
    sc.pp.scale(w,max_value=10);sc.tl.pca(w,n_comps=min(50,w.n_obs-1,w.n_vars-1),svd_solver='arpack',random_state=20260926)
    z=w.obsm['X_pca'];a.obsm['X_pca']=z
    if a.obs.study_technology.nunique()>1:
     h=np.asarray(hm.run_harmony(z,a.obs,'study_technology',random_state=20260926,ncores=4,verbose=False).Z_corr);h=h.T if h.shape!=z.shape else h
    else:h=z
    assert h.shape==z.shape and np.isfinite(h).all()
    a.obsm['X_harmony']=h;sc.pp.neighbors(a,use_rep='X_harmony',n_neighbors=min(15,a.n_obs-1),random_state=20260926)
    for key,res in [('cluster',1),('cluster_res05',.5)]:sc.tl.leiden(a,resolution=res,flavor='igraph',n_iterations=2,directed=False,key_added=key,random_state=20260926)
    sc.tl.umap(a,random_state=20260926);assert np.isfinite(a.obsm['X_umap']).all()
    cross=pd.crosstab(a.obs.original_cluster,a.obs.cluster);cross.to_csv(O/'tables'/f'{tag}_original_cluster_cross.tsv',sep='\t')
    cc=pd.crosstab(a.obs.frozen_class_id,a.obs.cluster);cc.to_csv(O/'tables'/f'{tag}_frozen_class_cross.tsv',sep='\t')
    for k,row in cc.iterrows():
     for cl,n in row.items():correspond.append(dict(layer=layer,scenario=scenario,frozen_class=k,new_cluster=cl,cells=int(n),fraction_of_class=float(n/row.sum())))
    for cl,sub in a.obs.groupby('cluster',observed=True):
     origgenes=['ACTA2','TAGLN','MYH11','MYL9','TPM2','CNN1','RGS5','CSPG4','NOTCH3','DES','MCAM','PDGFRB','COL1A1','COL1A2','SPP1','IBSP','RUNX2','SP7','PDGFRA','DPT','PI16','SFRP2']
     ix=a.obs.index.get_indexer(sub.index)
     for g in origgenes:
      j=np.flatnonzero(a.var.symbol.eq(g).to_numpy())
      if not len(j):continue
      v=np.asarray(a.X[ix,:][:,j].sum(1)).ravel()
      for donor,dsub in sub.groupby('donor_record',observed=True):
       jj=sub.index.get_indexer(dsub.index);profiles.append(dict(layer=layer,scenario=scenario,cluster=cl,donor_record=donor,gene=g,mean_logCP10k=float(v[jj].mean()),cells=len(jj),source='corrected_unregressed_expression' if scenario.startswith('decontx') else 'raw_unregressed_expression'))
    fig,axs=plt.subplots(1,3,figsize=(16,4.8),layout='constrained')
    for ax,col in zip(axs,['cluster','sample','frozen_class_id']):
     codes=pd.Categorical(a.obs[col]);ax.scatter(a.obsm['X_umap'][:,0],a.obsm['X_umap'][:,1],c=codes.codes,s=3,cmap='turbo',rasterized=True);ax.set_title(col+' ('+str(len(codes.categories))+')');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
     if len(codes.categories)<=24:
      for cat in codes.categories:
       v=a.obsm['X_umap'][a.obs[col].astype(str).eq(str(cat))];ax.text(*np.median(v,axis=0),str(cat),fontsize=7)
    fig.suptitle(tag);fig.savefig(O/'figures'/f'{tag}_UMAP.png',dpi=130);plt.close(fig)
    a.obs.to_csv(O/'tables'/f'{tag}_labels.tsv.gz',sep='\t',index_label='cell_id')
    q={'layer':layer,'scenario':scenario,'cells':a.n_obs,'same_cell_ids_and_order':list(a.obs_names)==list(base.obs_names),'finite_embedding':True,'clusters':a.obs.cluster.nunique(),'baseline_clusters':base.obs.cluster.nunique(),'ARI':float(adjusted_rand_score(a.obs.original_cluster,a.obs.cluster)),'NMI':float(normalized_mutual_info_score(a.obs.original_cluster,a.obs.cluster)),'ARI_resolution05':float(adjusted_rand_score(base.obs.cluster_res05.astype(str),a.obs.cluster_res05)),'no_biological_cells_removed':True,'no_threshold_changed':True}
    assert q['same_cell_ids_and_order']
    clean(a);a.write_h5ad(O/'objects'/f'{tag}.h5ad',compression='gzip');writej(done,q);summaries.append(q);outcomes.append(dict(layer=layer,scenario=scenario,status='complete'))
    pd.DataFrame(summaries).to_csv(O/'tables/sensitivity_summary.tsv',sep='\t',index=False)
    del a,w;gc.collect()
   except Exception as e:
    traceback.print_exc();outcomes.append(dict(layer=layer,scenario=scenario,status='failed',error=str(e)))
   writej(O/'tables/sensitivity_outcomes.json',outcomes)
  del base;gc.collect()
 pd.DataFrame(profiles).to_csv(O/'tables/sensitivity_cluster_donor_marker_profiles.tsv.gz',sep='\t',index=False)
 pd.DataFrame(correspond).to_csv(O/'tables/sensitivity_frozen_class_dispersion.tsv',sep='\t',index=False)
 qa={'at':now(),'expected_models':20,'recorded':len(outcomes),'complete':sum(r['status'].startswith('complete') for r in outcomes),'failures':[r for r in outcomes if r['status']=='failed'],'all_input_cells_retained':True,'classification_not_overwritten':True}
 writej(O/'delivery/SENSITIVITY_QA.json',qa);rc=0 if qa['complete']==20 else 1
 (O/'sensitivity.exit_code').write_text(str(rc));(O/'SENSITIVITY_FINISHED').write_text(now());status('MODELS_FINISHED_REVIEW_REQUIRED',**qa)
except Exception as e:
 traceback.print_exc();(O/'sensitivity.exit_code').write_text('1');status('FAILED',error=str(e));raise
