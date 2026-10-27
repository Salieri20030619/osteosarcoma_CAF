from pathlib import Path
import os,sys,json,hashlib,datetime,gc,gzip,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='8'
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse,io
from sklearn.metrics import adjusted_rand_score
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase53_de_novo_reannotation_20261002'
dataset=sys.argv[1];R=O/'objects'/dataset
for d in [R,O/'tables'/dataset,O/'figures'/dataset]:d.mkdir(parents=True,exist_ok=True)
seed=20261002;np.random.seed(seed)
m=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);m=m[m.dataset.eq(dataset)]
cw=pd.read_csv(P/'phase43_full_transcriptome_discovery_20261001/tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False)
symbols=dict(zip(cw.feature,cw.canonical_symbol));quality=dict(zip(cw.feature,~cw.mapping_method.eq('HGNC_Ensembl_approved')))
parts=[];sets=[];ledger=[];libs=[];features_by_sample={}
for row in m.itertuples(index=False):
 f=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad'
 a=ad.read_h5ad(f);x=a.X.tocsr().astype(np.float32)
 obs=pd.DataFrame(index=a.obs_names.astype(str));obs['sample']=row.sample;obs['dataset']=dataset;obs['patient_id']=row.donor_record;obs['site']=row.site;obs['technology']=row.technology
 lib=np.asarray(x.sum(1)).ravel().astype(float);obs['full_library']=lib
 hold=np.array([symbols.get(str(g),'') in ['COL1A1','LOX','MMP14'] for g in a.var_names])
 obs['analysis_library']=lib-np.asarray(x[:,hold].astype(float).sum(1)).ravel()
 obs['nonpositive_analysis_library']=obs.analysis_library.le(0).to_numpy()
 parts.append(ad.AnnData(x,obs=obs,var=pd.DataFrame(index=a.var_names.astype(str))))
 sets.append(set(a.var_names.astype(str)));features_by_sample[row.sample]=a.n_vars
 ledger.append({'sample':row.sample,'dataset':dataset,'cells':a.n_obs,'features':a.n_vars,'input_SHA256':hashlib.sha256(f.read_bytes()).hexdigest()})
 del a,x;gc.collect()
common=sorted(set.intersection(*sets));allfeatures=sorted(set.union(*sets));assert len(common)>1000
# Same-source common assay features avoid treating unmeasured genes as biological zero.
full=ad.concat([a[:,common] for a in parts],join='inner',merge='same',index_unique=None);del parts;gc.collect()
assert full.obs_names.is_unique and len(full)==sum(x['cells'] for x in ledger)
raw=full.X.tocsr();obs=full.obs.copy();var=list(full.var_names)
canon=np.array([symbols.get(g,'') for g in var],dtype=object)
hold=np.isin(canon,['COL1A1','LOX','MMP14'])
den=np.where(obs.analysis_library.to_numpy()>0,obs.analysis_library.to_numpy(),1.)
norm=raw.astype(np.float64).multiply((1e4/den)[:,None]).tocsr();norm.data=np.log1p(norm.data)
model=ad.AnnData(norm[:,~hold].astype(np.float32),obs=obs.copy(),var=pd.DataFrame(index=np.array(var)[~hold]))
sc.pp.highly_variable_genes(model,flavor='seurat',n_top_genes=min(3000,model.n_vars),batch_key='sample',subset=False)
pd.DataFrame({'feature':model.var_names,'canonical_symbol':[symbols.get(g,'') for g in model.var_names],'highly_variable':model.var.highly_variable.to_numpy()}).to_csv(O/'tables'/dataset/'structure_features.tsv',sep='\t',index=False)
model=model[:,model.var.highly_variable].copy();assert not set(symbols.get(g,'') for g in model.var_names)&set(['COL1A1','LOX','MMP14'])
sc.pp.scale(model,zero_center=False,max_value=None)
sc.tl.pca(model,n_comps=min(50,model.n_vars-1,model.n_obs-1),zero_center=True,svd_solver='arpack',random_state=seed)
if model.obs['sample'].nunique()>1:
 h=hm.run_harmony(model.obsm['X_pca'],model.obs[['sample']],['sample'],max_iter_harmony=20,random_state=seed,ncores=8,verbose=False)
 emb=np.asarray(h.Z_corr);emb=emb if emb.shape[0]==model.n_obs else emb.T
 assert emb.shape==model.obsm['X_pca'].shape
 model.obsm['X_integrated']=emb.astype(np.float32)
else:model.obsm['X_integrated']=model.obsm['X_pca'].copy()
sc.pp.neighbors(model,n_neighbors=15,use_rep='X_integrated',random_state=seed)
for res,key in [(.5,'r05'),(1.,'r10'),(2.,'r20')]:
 sc.tl.leiden(model,resolution=res,key_added=key,flavor='igraph',directed=False,n_iterations=2,random_state=seed)
 model.obs[key]=model.obs[key].astype(str)
sc.tl.umap(model,random_state=seed,maxiter=200)
for key in ['r05','r10','r20']:obs[key]=model.obs[key].to_numpy()
obs['UMAP1']=model.obsm['X_umap'][:,0];obs['UMAP2']=model.obsm['X_umap'][:,1]
obs.to_csv(O/'tables'/dataset/'cell_structure.tsv.gz',sep='\t',index_label='cell_id')
np.savez_compressed(R/'embeddings.npz',cell_id=obs.index.to_numpy(str),PCA=model.obsm['X_pca'],integrated=model.obsm['X_integrated'],UMAP=model.obsm['X_umap'])
sparse.save_npz(R/'connectivities.npz',model.obsp['connectivities'])
# Export global clusters at all resolutions and patient clusters at the primary resolution.
profiles=[];meta=[]
def add_profiles(labels,kind,resolution,patient_labels=None):
 unique=sorted(set(labels));index={g:j for j,g in enumerate(unique)};codes=np.array([index[g] for g in labels])
 gmat=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(unique),len(codes)))
 sums=(gmat@raw).tocsr();l=np.asarray(gmat@obs.analysis_library.to_numpy()).ravel()
 for j,g in enumerate(unique):
  selected=codes==j;prof_id=dataset+'__'+kind+'__'+resolution+'__'+str(j)
  profiles.append(sums[j]);meta.append({'profile_id':prof_id,'kind':kind,'resolution':resolution,'cluster':str(g) if kind=='global' else str(g).split('|||')[-1],'patient_id':'' if kind=='global' else str(g).split('|||')[0],'cells':int(selected.sum()),'patients':obs.loc[selected,'patient_id'].nunique(),'samples':obs.loc[selected,'sample'].nunique(),'analysis_library':float(l[j]),'max_patient_fraction':float(obs.loc[selected,'patient_id'].value_counts(normalize=True).max())})
for key in ['r05','r10','r20']:add_profiles(obs[key].astype(str).to_numpy(),'global',key)
add_profiles((obs.patient_id.astype(str)+'|||'+obs.r10.astype(str)).to_numpy(),'patient','r10')
pb=sparse.vstack(profiles,format='csr');pm=pd.DataFrame(meta)
# Canonical feature choice is made on the assay universe, independently of outcome.
choices={}
for j,f in enumerate(var):
 g=symbols.get(f,'');k=(quality.get(f,1),f)
 if g and g not in ['COL1A1','LOX','MMP14'] and (g not in choices or k<choices[g][0]):choices[g]=(k,j)
ann_genes=sorted(choices);ann_idx=[choices[g][1] for g in ann_genes]
ann=pb[:,ann_idx].toarray().T.astype(np.float64);ld=np.where(pm.analysis_library.to_numpy()>0,pm.analysis_library.to_numpy(),1)
ann=np.log1p(ann/ld[None,:]*1e4)
import h5py
with h5py.File(R/'reference_input.h5','w') as f:
 f.create_dataset('logcounts',data=ann,compression='gzip')
 f.create_dataset('genes',data=np.asarray(ann_genes,dtype='S'))
 f.create_dataset('profiles',data=np.asarray(pm.profile_id,dtype='S'))
pm.to_csv(O/'tables'/dataset/'reference_profile_metadata.tsv',sep='\t',index=False)
# Export full source common-gene expression for all primary clusters, independent of proposed labels.
globalpm=pm[(pm.kind=='global')&(pm.resolution=='r10')]
ri=globalpm.index.to_numpy();cl_ids=globalpm.cluster.tolist();cl_pb=pb[ri]
io.mmwrite(str(R/'primary_cluster_counts.mtx'),cl_pb.T); 
with open(R/'primary_cluster_counts.mtx','rb') as f,gzip.open(R/'primary_cluster_counts.mtx.gz','wb') as z:
 import shutil;shutil.copyfileobj(f,z)
(R/'primary_cluster_counts.mtx').unlink()
pd.DataFrame({'feature':var,'symbol':canon}).to_csv(O/'tables'/dataset/'genes.tsv',sep='\t',index=False)
globalpm.to_csv(O/'tables'/dataset/'primary_clusters.tsv',sep='\t',index=False)
# Descriptive marker evidence: same-patient, same-site cluster-versus-rest, not cell-level significance tests.
loc=obs.patient_id.astype(str)+'|||'+obs.site.astype(str);locs=sorted(set(loc));cr=[]
clusterindex={k:i for i,k in enumerate(cl_ids)}
locinfo={}
for lk in locs:
 base=(loc==lk).to_numpy();locinfo[lk]=(base,np.asarray(raw[base].sum(0,dtype=np.float64)).ravel(),float(obs.loc[base,'full_library'].sum()))
total_detect=np.asarray((raw>0).sum(0)).ravel();total_mean=np.asarray(norm.sum(0)).ravel()/raw.shape[0]
cldesc=[]
for cl in cl_ids:
 mask=obs.r10.eq(cl).to_numpy();n=int(mask.sum());counts=np.asarray(raw[mask].sum(0)).ravel();det=np.asarray((raw[mask]>0).sum(0)).ravel()/n
 mn=np.asarray(norm[mask].sum(0)).ravel()/n;restmn=(total_mean*len(obs)-mn*n)/(len(obs)-n) if n<len(obs) else np.full(len(var),np.nan)
 patientsite=[]
 for lk in locs:
  base,total_counts,total_lib=locinfo[lk];left=base&mask;right=base&~mask
  if not left.any() or not right.any():continue
  ll=obs.loc[left,'full_library'].sum();rr=total_lib-ll
  if ll<=0 or rr<=0:continue
  a=np.asarray(raw[left].sum(0,dtype=np.float64)).ravel();b=total_counts-a
  d=np.log2(a/ll*1e4+.1)-np.log2(b/rr*1e4+.1)
  patientsite.append((lk.split('|||')[0],d))
 vectors=[]
 for patient in sorted(set(x[0] for x in patientsite)):vectors.append(np.mean([x[1] for x in patientsite if x[0]==patient],axis=0))
 if vectors:
  vv=np.stack(vectors);delta=vv.mean(0);positive=(vv>0).sum(0);negative=(vv<0).sum(0);N=len(vectors)
 else:delta=np.full(len(var),np.nan);positive=np.zeros(len(var),int);negative=np.zeros(len(var),int);N=0
 table=pd.DataFrame({'dataset':dataset,'cluster':cl,'feature':var,'symbol':canon,'cells':n,'N_patients_matched':N,'mean_log_expression':mn,'rest_mean_log_expression':restmn,'detection_fraction':det,'patient_mean_log2_difference':delta,'positive_patients':positive,'negative_patients':negative,'evidence_type':'descriptive_annotation_no_P'})
 table.to_csv(O/'tables'/dataset/f'cluster_{cl}_full_marker_evidence.tsv.gz',sep='\t',index=False)
 # Show observed positive expression and rank effect; this is descriptive ranking, not a hypothesis screen.
 top=table[table.symbol.ne('')&table.detection_fraction.gt(0)].sort_values(['patient_mean_log2_difference','detection_fraction'],ascending=[False,False]).drop_duplicates('symbol').head(40)
 top.to_csv(O/'tables'/dataset/f'cluster_{cl}_top_markers.tsv',sep='\t',index=False)
 cldesc.append({'dataset':dataset,'cluster':cl,'cells':n,'patients':obs.loc[mask,'patient_id'].nunique(),'N_patients_matched':N,'max_patient_fraction':float(obs.loc[mask,'patient_id'].value_counts(normalize=True).max()),'top_markers':';'.join(top.symbol.head(20))})
 print('MARKERS_COMPLETE',dataset,cl,n,N,flush=True)
pd.DataFrame(cldesc).to_csv(O/'tables'/dataset/'cluster_evidence_summary.tsv',sep='\t',index=False)
pd.DataFrame(ledger).to_csv(O/'tables'/dataset/'input_manifest.tsv',sep='\t',index=False)
qa={'dataset':dataset,'cells':len(obs),'samples':len(m),'patients':obs.patient_id.nunique(),'source_features_union':len(allfeatures),'source_features_common':len(common),'heldout_features':int(hold.sum()),'structure_HVG':model.n_vars,'clusters':{k:int(obs[k].nunique()) for k in ['r05','r10','r20']},'resolution_ARI_05_10':float(adjusted_rand_score(obs.r05,obs.r10)),'resolution_ARI_10_20':float(adjusted_rand_score(obs.r10,obs.r20)),'reference_profiles':len(pm),'reference_genes':len(ann_genes),'no_old_labels_used':True,'no_added_cell_filter':True,'nonpositive_analysis_library':int(obs.nonpositive_analysis_library.sum())}
(O/'tables'/dataset/'STRUCTURE_QA.json').write_text(json.dumps(qa,indent=2));(R/'STRUCTURE_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
