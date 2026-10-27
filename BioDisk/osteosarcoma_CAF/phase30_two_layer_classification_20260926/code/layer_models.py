from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,fcntl,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase30_two_layer_classification_20260926'
for d in ['code','resources','tables','objects','figures','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):
 t=f.with_name(f.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS.json',{'status':s,'at':now(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def clean(a):
 for frame in [a.obs,a.var]:
  for c in frame:
   if frame[c].dtype==object or isinstance(frame[c].dtype,pd.StringDtype):frame[c]=frame[c].fillna('').astype(str).astype(object)
  frame.index=pd.Index(frame.index.astype(str).to_numpy(dtype=object))
def model(a,layer):
 status('FITTING_LAYER',layer=layer,cells=a.n_obs)
 a.layers['counts']=a.X.copy()
 a.X=a.X.astype(np.float32).multiply((1e4/a.obs.full_library_counts.to_numpy())[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
 sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='study_technology')
 w=a[:,a.var.highly_variable].copy();sc.pp.scale(w,max_value=10)
 sc.tl.pca(w,n_comps=min(50,w.n_obs-1,w.n_vars-1),svd_solver='arpack',random_state=20260926)
 z=w.obsm['X_pca'];a.obsm['X_pca']=z
 if a.obs.study_technology.nunique()>1:
  h=hm.run_harmony(z,a.obs,'study_technology',random_state=20260926,ncores=4,verbose=False)
  h=np.asarray(h.Z_corr);h=h.T if h.shape!=z.shape else h;assert h.shape==z.shape and np.isfinite(h).all()
 else:h=z
 a.obsm['X_harmony']=h;sc.pp.neighbors(a,use_rep='X_harmony',n_neighbors=min(15,a.n_obs-1),random_state=20260926)
 sc.tl.leiden(a,resolution=1,flavor='igraph',n_iterations=2,directed=False,key_added='cluster',random_state=20260926)
 sc.tl.leiden(a,resolution=.5,flavor='igraph',n_iterations=2,directed=False,key_added='cluster_res05',random_state=20260926)
 sc.tl.umap(a,random_state=20260926)
 del w;gc.collect();return a
def profile(a,layer):
 rows=[];markers=[];programs=[];cols=[c for c in a.obs if c.startswith(('lineage_','program_'))]
 for cl,g in a.obs.groupby('cluster',observed=True):
  ix=a.obs.index.get_indexer(g.index);mask=np.zeros(a.n_obs,bool);mask[ix]=True
  dn=g.groupby('donor_record',observed=True).size()
  row={'layer':layer,'cluster':str(cl),'cells':len(g),'samples':g['sample'].nunique(),'donor_records':len(dn),'datasets':g.dataset.nunique(),'largest_donor_fraction':float(dn.max()/len(g)),'largest_donor':dn.idxmax(),'supported_cells':int(g.in_supported.sum()),'dual_reference_cells':int(g.dual_reference.sum()),'identity_unresolved_cells':int((~g.in_supported).sum())}
  for c in cols:row[c]=float(g.groupby('donor_record',observed=True)[c].mean().mean())
  rows.append(row)
  for (ds,dn,site),q in g.groupby(['dataset','donor_record','site'],observed=True):
   programs.append(dict(layer=layer,cluster=str(cl),dataset=ds,donor_record=dn,site=site,cells=len(q),**{c:float(q[c].mean()) for c in cols}))
  mean=np.asarray(a.X[mask].mean(0)).ravel();other=np.asarray(a.X[~mask].mean(0)).ravel()
  frac=a.layers['counts'][mask].getnnz(0)/mask.sum();fraco=a.layers['counts'][~mask].getnnz(0)/(~mask).sum()
  for j in np.argsort(mean-other)[::-1][:80]:
   markers.append(dict(layer=layer,cluster=str(cl),ensembl=a.var_names[j],symbol=a.var.symbol.iloc[j],delta_mean_logCP10k=float(mean[j]-other[j]),mean_logCP10k=float(mean[j]),detection_fraction=float(frac[j]),rest_detection_fraction=float(fraco[j]),evidence='descriptive_marker_ranking_not_inferential'))
 save(pd.DataFrame(rows),layer+'_cluster_profiles.tsv');save(pd.DataFrame(markers),layer+'_descriptive_markers.tsv');save(pd.DataFrame(programs),layer+'_donor_site_programs.tsv')
 return pd.DataFrame(rows)
def plot(a,layer,pr):
 fig,axs=plt.subplots(1,3,figsize=(17,5),layout='constrained')
 for ax,col,title in zip(axs,['cluster','dataset','in_supported'],['Expression clusters','Dataset','Supported candidate membership']):
  for name,g in a.obs.groupby(col,observed=True):
   ix=a.obs.index.get_indexer(g.index);ax.scatter(a.obsm['X_umap'][ix,0],a.obsm['X_umap'][ix,1],s=3,label=str(name),rasterized=True)
  ax.set_title(title);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(fontsize=6,markerscale=2,frameon=False)
 fig.suptitle(layer+' candidate layer');fig.savefig(O/'figures'/f'{layer}_UMAP.png',dpi=150);plt.close(fig)
 cols=[c for c in pr if c.startswith(('program_','lineage_'))];v=pr[cols].to_numpy();z=(v-v.mean(0))/np.maximum(v.std(0),1e-8)
 fig,ax=plt.subplots(figsize=(14,max(4,len(pr)*.35)),layout='constrained')
 im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(len(pr)),pr.cluster);ax.set_xticks(range(len(cols)),[c.replace('lineage_','L:').replace('program_','P:') for c in cols],rotation=70,ha='right');fig.colorbar(im,ax=ax,label='Standardized donor-mean program score')
 ax.set_title(layer+' marker review');fig.savefig(O/'figures'/f'{layer}_programs.png',dpi=150);plt.close(fig)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'LAYER_MODELS_COMPLETE').exists()
rc=1
try:
 src=P/'phase28_CAF_clustering_20260926/objects/CAF_integrated.h5ad'
 a=ad.read_h5ad(src);assert a.n_obs==11379
 ann=pd.read_csv(P/'phase28_CAF_clustering_20260926/tables/CAF_cluster_annotation_reviewed.tsv',sep='\t').set_index('cluster');ann.index=ann.index.astype(str)
 current=pd.read_csv(P/'phase29_myCAF_followup_20260926/tables/CAF_reviewed_labels.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
 supported=current.index[current.identity_review.eq('retain_CAF_candidate')]
 a.obs['source_cluster28']=a.obs.cluster.astype(str)
 clear={'1','5','6','9','11','12','15','19','24','26'}
 a.obs['in_extended']=~a.obs.source_cluster28.isin(clear)
 a.obs['in_supported']=a.obs_names.isin(supported)
 a.obs['dual_reference']=a.obs.in_supported&a.obs.CNV_evidence.eq('concordant_no_alteration')
 a.obs['previous_lineage_review']=a.obs.source_cluster28.map(ann.CAF_type)
 a.obs['candidate_identity_status']=np.where(a.obs.in_supported,'supported_candidate',np.where(a.obs.in_extended,'identity_unresolved_exploratory','competing_lineage_background'))
 a.obs['previous_supported_type']=current.CAF_type.reindex(a.obs_names).fillna('not_in_previous_refined_pool')
 a.obs['clinical_group']=a.obs['sample'].map(pd.read_csv(P/'phase27_identity_review_20260926/tables/clinical_mapping_v3.tsv',sep='\t').fillna('').set_index('sample').clinical_group)
 assert int(a.obs.in_extended.sum())==8174 and int(a.obs.in_supported.sum())==1224 and int(a.obs.dual_reference.sum())==114
 assert (a.obs.in_extended|~a.obs.in_supported).all()
 a.obs.to_csv(O/'tables/layer_membership_all11379.tsv.gz',sep='\t',index_label='cell_id')
 protocol={'created':now(),'user_scope':'Main research in supported1224 and expanded8174; dual114 is evidence reference only, no separate clustering','layers_nested':True,'expanded_unresolved_cells':6950,'competing_lineage_background':3205,'expanded_identity':'Candidate expression space including osteogenic/mitochondrial/cycling identity-uncertain cells; never equated with confirmed CAF','no_new_cell_entry_thresholds':True,'model':'Independent layer refits: shared genes, full-library logCP10k, HVG2000, PCA50, Harmony study_technology defaults, neighbors15, Leiden1 primary /0.5 boundary sensitivity; fixed seed20260926','classification':'Assign biological names from marker and program profiles BEFORE inferential statistics; preserve technical/identity uncertainty flags','statistics_planned':'Patient-level within-layer class vs other candidates, stratified by cohort and matched sampling site; nominal P<0.05 may nominate exploratory candidates; complete per-family P/FDR retained. No parameter changes to seek significance. No cell-level biological inference.','source_frozen':str(src),'active_S525_recovery':'Independent; not incorporated into frozen layers without a separate review','old_results_preserved':True}
 writej(O/'PROTOCOL.json',protocol);writej(O/'resources/SOURCE_SNAPSHOT.json',{'path':str(src),'bytes':src.stat().st_size,'mtime_ns':src.stat().st_mtime_ns,'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
 all_summary=[];models={}
 for layer,mask in [('supported',a.obs.in_supported),('expanded',a.obs.in_extended)]:
  b=a[mask].copy();b.X=b.layers['counts'].copy()
  b=model(b,layer);pr=profile(b,layer);plot(b,layer,pr);clean(b)
  b.write_h5ad(O/'objects'/f'{layer}_integrated.h5ad',compression='gzip')
  b.obs.to_csv(O/'tables'/f'{layer}_cell_labels.tsv.gz',sep='\t',index_label='cell_id')
  models[layer]=b
  all_summary.append({'layer':layer,'cells':b.n_obs,'samples':b.obs['sample'].nunique(),'donor_records_not_pooled_independent_N':b.obs.donor_record.nunique(),'datasets':b.obs.dataset.nunique(),'clusters':b.obs.cluster.nunique(),'dual_reference_cells':int(b.obs.dual_reference.sum()),'finite_embedding':bool(np.isfinite(b.obsm['X_umap']).all())})
 s=models['supported'];e=models['expanded'];cross=pd.crosstab(s.obs.cluster,e.obs.loc[s.obs_names,'cluster'])
 cross.index.name='supported_cluster';cross.columns.name='expanded_cluster';cross.to_csv(O/'tables/cluster_correspondence_counts.tsv',sep='\t')
 cross.div(cross.sum(axis=1),axis=0).to_csv(O/'tables/cluster_correspondence_row_fraction.tsv',sep='\t')
 qa={'at':now(),'layers':all_summary,'shared_supported_cells':s.n_obs,'ARI_on_shared_cells':float(adjusted_rand_score(s.obs.cluster,e.obs.loc[s.obs_names,'cluster'])),'114_reference_separately_clustered':False,'source_counts_not_changed':True,'next':'Review profiles, assign new class labels, then patient-level marker and clinical exploration with nominal significance explicitly labelled'}
 writej(O/'delivery/LAYER_QA.json',qa);save(pd.DataFrame(all_summary),'layer_summary.tsv')
 (O/'LAYER_MODELS_COMPLETE').write_text(now()+chr(10));status('TWO_LAYER_MODELS_COMPLETE_CLASSIFICATION_REVIEW',**qa);rc=0
except Exception as ex:status('FAILED',error=str(ex));traceback.print_exc()
finally:(O/'models.exit_code').write_text(str(rc)+chr(10))
raise SystemExit(rc)
