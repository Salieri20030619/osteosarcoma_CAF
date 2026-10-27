from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921';SEED=20260921
raise SystemExit('SUPERSEDED: customer requires CNV and lineage identity screening before CAF subtype modeling; all-cell model is not the CAF main analysis.')
assert (O/'INPUTS_COMPLETE').exists()
def status(s,**kw):
 (O/'MODEL_STATUS.json').write_text(json.dumps({'status':s,'updated':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw},indent=2));print(s,kw,flush=True)
def truth(s):return s.fillna(False).astype(str).str.lower().isin(['true','1','1.0'])
def clean(a):
 for frame in [a.obs,a.var]:
  for c in frame:
   if frame[c].dtype==object or isinstance(frame[c].dtype,pd.StringDtype):frame[c]=frame[c].fillna('missing').astype(str).astype(object)
  frame.index=pd.Index(frame.index.to_numpy(dtype=object),dtype=object)
def save(a,path):clean(a);tmp=path.with_suffix('.tmp.h5ad');a.write_h5ad(tmp,compression='gzip');tmp.replace(path)
done=O/'objects/all_entry_integrated.h5ad'
if done.exists() and (O/'MODEL_COMPLETE').exists():status('MODEL_CACHED');raise SystemExit(0)
checkpoint=O/'objects/03_neighbors.h5ad'
if checkpoint.exists():a=ad.read_h5ad(checkpoint)
else:
 checkpoint2=O/'objects/02_harmony.h5ad'
 if checkpoint2.exists():a=ad.read_h5ad(checkpoint2)
 else:
  checkpoint1=O/'objects/01_pca.h5ad'
  if checkpoint1.exists():a=ad.read_h5ad(checkpoint1)
  else:
   hvg=(O/'resources/HVG2000.txt').read_text().splitlines();meta=pd.read_csv(O/'tables/input_recovery.tsv',sep='\t');parts=[];obs=[];symbols={}
   status('LOADING_ALL_ENTRY_HVG',samples=len(meta))
   for sm in meta['sample']:
    b=ad.read_h5ad(O/'inputs'/f'{sm}.h5ad',backed='r');inds=pd.Series(np.arange(b.n_vars),index=b.var.ensembl_id.astype(str));assert inds.index[inds.index.str.startswith('ENSG')].is_unique
    ix=inds.loc[hvg].to_numpy();mask=truth(b.obs.embeddable).to_numpy();x=b.X[:,np.sort(ix)].tocsr();x=x[:,np.argsort(np.argsort(ix))][mask].astype(np.float32)
    den=b.obs.n_counts.to_numpy()[mask];x=x.multiply((1e4/den)[:,None]).tocsr();x.data=np.log1p(x.data)
    ob=b.obs.loc[mask].copy();ob['known_patient']=~ob.patient_id.fillna('').astype(str).str.lower().isin(['','nan','na','unknown','unconfirmed','not_established','missing'])
    ob['donor_record']=np.where(ob.known_patient,ob.dataset.astype(str)+'|'+ob.patient_id.astype(str),ob.dataset.astype(str)+'|'+ob['sample'].astype(str)+'|UNKNOWN_PATIENT')
    parts.append(x);obs.append(ob);symbols.update(dict(zip(b.var.ensembl_id,b.var.symbol)));b.file.close()
   a=ad.AnnData(sparse.vstack(parts).tocsr(),obs=pd.concat(obs),var=pd.DataFrame({'symbol':[symbols[g] for g in hvg]},index=hvg));assert a.obs_names.is_unique
   del parts,obs;gc.collect();status('PCA_RUNNING',cells=a.n_obs,genes=a.n_vars)
   # Scaling and dimensionality choices affect the model, never membership.
   w=a.copy();sc.pp.scale(w,max_value=10);sc.tl.pca(w,n_comps=30,svd_solver='arpack',random_state=SEED);a.obsm['X_pca']=w.obsm['X_pca'];a.varm['PCs']=w.varm['PCs'];a.uns['pca']=w.uns['pca'];del w;gc.collect();save(a,checkpoint1)
  status('HARMONY_RUNNING',cells=a.n_obs)
  z=a.obsm['X_pca'];h=hm.run_harmony(z,a.obs,'study_technology',theta=2,random_state=SEED,ncores=4,max_iter_harmony=30,verbose=False);zz=np.asarray(h.Z_corr);zz=zz.T if zz.shape!=z.shape else zz;assert zz.shape==z.shape and np.isfinite(zz).all();a.obsm['X_harmony']=zz;save(a,checkpoint2)
 status('NEIGHBORS_AND_CLUSTERS_RUNNING',cells=a.n_obs)
 sc.pp.neighbors(a,use_rep='X_harmony',n_neighbors=15,random_state=SEED)
 sc.tl.leiden(a,resolution=.5,key_added='cluster',flavor='igraph',n_iterations=2,directed=False,random_state=SEED);save(a,checkpoint)
status('UMAP_RUNNING',cells=a.n_obs,clusters=a.obs.cluster.nunique());sc.tl.umap(a,random_state=SEED)
prior=pd.read_csv(P/'phase17_all_cohort_expanded_20260921/tables/supported_cell_labels.tsv.gz',sep='\t',dtype={'cluster':str}).set_index('cell_id')
a.obs['prior_supported_cluster']=prior.cluster.reindex(a.obs_names).fillna('not_in_old_supported').to_numpy()
save(a,done);a.obs.to_csv(O/'tables/new_all_entry_labels.tsv.gz',sep='\t',index_label='cell_id')
scores=[c for c in a.obs if c.startswith('score_')];donors=a.obs.groupby(['cluster','donor_record'],observed=True)[scores].mean();prof=donors.groupby('cluster',observed=True).mean();prof.to_csv(O/'tables/cluster_donor_balanced_programs.tsv',sep='\t')
rows=[]
for cl,g in a.obs.groupby('cluster',observed=True):
  d=g.groupby('donor_record',observed=True).size();rows.append({'cluster':cl,'cells':len(g),'samples':g['sample'].nunique(),'known_donor_records':g.loc[truth(g.known_patient),'donor_record'].nunique(),'display_records':g.donor_record.nunique(),'datasets':g.dataset.nunique(),'max_record_fraction':d.max()/len(g),'median_genes':g.n_genes.median(),'median_UMI':g.n_counts.median(),'median_MT':g.pct_mt.median(),'restored_cells':int((~truth(g.historical_QC_keep)).sum()),'old_strict':int(truth(g.historical_strict).sum()),'old_supported':int(truth(g.historical_supported).sum()),'old_C7':int(g.prior_supported_cluster.eq('7').sum()),'old_C9':int(g.prior_supported_cluster.eq('9').sum())})
pd.DataFrame(rows).to_csv(O/'tables/new_cluster_summary.tsv',sep='\t',index=False)
pd.crosstab(a.obs.cluster,a.obs.prior_supported_cluster).to_csv(O/'tables/old_supported_new_cluster_overlap.tsv',sep='\t')
if 'official_consensus_celltype_annotation' in a.obs:pd.crosstab(a.obs.cluster,a.obs.official_consensus_celltype_annotation).to_csv(O/'tables/official_annotation_cluster_overlap.tsv',sep='\t')
xy=a.obsm['X_umap'];rng=np.random.default_rng(SEED);ix=np.arange(a.n_obs);rng.shuffle(ix)
fig,axs=plt.subplots(1,3,figsize=(18,5))
for ax,c in zip(axs,['cluster','dataset','historical_QC_keep']):
 cat=a.obs[c].astype('category');ax.scatter(xy[ix,0],xy[ix,1],c=cat.cat.codes.to_numpy()[ix],s=.4,cmap='tab20',rasterized=True);ax.set(title=c,xticks=[],yticks=[])
 if c=='cluster':
  for cl in cat.cat.categories:
   mid=np.median(xy[a.obs.cluster.astype(str).eq(str(cl))],axis=0);ax.text(*mid,str(cl),fontsize=8,ha='center',bbox=dict(fc='white',alpha=.7,pad=.5,ec='none'))
fig.tight_layout();fig.savefig(O/'figures/all_entry_UMAP.png',dpi=180);plt.close(fig)
z=(prof-prof.mean())/prof.std().replace(0,1);fig,ax=plt.subplots(figsize=(12,max(6,len(z)*.24)));im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set(yticks=np.arange(len(z)),yticklabels=z.index,xticks=np.arange(len(z.columns)),xticklabels=[s.replace('score_','') for s in z.columns]);plt.setp(ax.get_xticklabels(),rotation=65,ha='right');fig.colorbar(im,ax=ax,label='Cluster program z score; donor-balanced');fig.tight_layout();fig.savefig(O/'figures/all_cluster_programs.png',dpi=180);plt.close(fig)
(O/'MODEL_COMPLETE').write_text('ALL_ENTRY_MODEL_COMPLETE\n');status('COMPLETE',cells=a.n_obs,clusters=a.obs.cluster.nunique())
