from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ['TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD']='1'
import json,datetime,traceback,gc,hashlib
import numpy as np,pandas as pd,anndata as ad,scvi,torch,h5py
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer
from lightning.pytorch import Callback
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase32_lung_sensitivity_20260926';M=P/'phase30_two_layer_classification_20260926';B=P/'phase26_client_reanalysis_20260922';A=P/'phase25_complete_nonPDX_assets_20260921/references/HLCA'
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):
 t=f.with_name(f.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'HLCA_STATUS.json',dict(status=s,at=now(),**kw));print(s,kw,flush=True)
def clean(a):
 for fr in [a.obs,a.var]:
  for c in fr:
   if fr[c].dtype==object or isinstance(fr[c].dtype,pd.StringDtype):fr[c]=fr[c].fillna('').astype(str).astype(object)
  fr.index=pd.Index(fr.index.astype(str).to_numpy(dtype=object))
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):status('QUERY_TRAINING',epoch=int(trainer.current_epoch)+1,maximum_epochs=500)
protocol={'at':now(),'method':'HLCA pretrained scANVI/scArches + distance-weighted 50-nearest-reference label transfer','reference_model':'original HLCA legacy0.8.1 saved model converted with official SCANVI.convert_legacy_save','reference_embedding':'HLCA_full_v1.1_emb.h5ad X is 30D latent embedding, core cells only; never used as expression counts','query':'8174 frozen expanded candidates; supported1224 reported as subset, not independent validation','input':'original integer counts in full measured source feature space projected to ordered reference2000; missing source genes padded to zero with complete coverage records','training':{'max_epochs':500,'early_stopping':True,'early_stopping_monitor':'elbo_train','early_stopping_patience':10,'early_stopping_min_delta':.001,'weight_decay':0.,'seed':20260926,'accelerator':'cpu','threads':4},'knn':{'n_neighbors':50,'metric':'euclidean','algorithm':'brute','weights':'official scArches distance-to-affinity formula','query_chunk_size_resource_only':256},'confidence_filter':None,'uncertainty_retained_continuously':True,'no_cells_dropped':True,'interpretation':'reference phenotype similarity, not evidence of normal origin, confirmed CAF identity or tumor specificity','official_sources':['https://docs.scvi-tools.org/en/stable/tutorials/notebooks/hub/query_hlca_knn.html','https://raw.githubusercontent.com/theislab/scarches/master/scarches/utils/knn.py']}
writej(O/'resources/HLCA_PROTOCOL.json',protocol)
try:
 scvi.settings.seed=20260926;torch.set_num_threads(4)
 modeldir=O/'resources/HLCA_converted';saved=torch.load(modeldir/'model.pt',map_location='cpu',weights_only=False);genes=pd.Index(saved['var_names'].astype(str));assert len(genes)==2000 and genes.is_unique
 ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id');assert len(ev)==8174
 su=read(M/'tables/supported_classified_cells.tsv.gz').set_index('cell_id');assert len(su)==1224
 queryfile=O/'objects/HLCA_query_raw2000.h5ad'
 if not queryfile.exists():
  mats=[];ids=[];coverage=[]
  for sm,sub in ev.groupby('sample',sort=False):
   status('ASSEMBLING_QUERY',sample=sm,cells_assembled=len(ids))
   a=ad.read_h5ad(B/'inputs'/sm/'full_called_counts.h5ad',backed='r');gi=a.var_names.get_indexer(genes);ok=gi>=0
   ci=a.obs_names.get_indexer(sub.index);assert (ci>=0).all()
   x=a[ci,:].X[:,gi[ok]].tocoo()
   mat=sparse.csr_matrix((x.data,(x.row,np.flatnonzero(ok)[x.col])),shape=(len(sub),len(genes)))
   mats.append(mat.astype(np.float32));ids.extend(sub.index)
   coverage.append(dict(sample=sm,cells=len(sub),reference_genes=2000,measured=int(ok.sum()),missing=int((~ok).sum()),missing_gene_ids=';'.join(genes[~ok]),no_coverage_gate=True));a.file.close()
  q=ad.AnnData(sparse.vstack(mats,format='csr'),obs=ev.loc[ids].copy(),var=pd.DataFrame(index=genes));q=q[ev.index].copy()
  q.obs['source_dataset']=q.obs.dataset.astype(str);q.obs['dataset']='CAF_QUERY__'+q.obs.study_technology.astype(str);q.obs['scanvi_label']='unlabeled'
  assert np.isfinite(q.X.data).all() and (q.X.data>=0).all() and np.equal(q.X.data,np.floor(q.X.data)).all()
  clean(q);q.write_h5ad(queryfile,compression='gzip');pd.DataFrame(coverage).to_csv(O/'tables/HLCA_query_gene_coverage.tsv',sep='\t',index=False)
 else:q=ad.read_h5ad(queryfile)
 status('LOADING_QUERY_MODEL',query_cells=q.n_obs,reference_genes=q.n_vars)
 scvi.model.SCANVI.prepare_query_anndata(q,str(modeldir));model=scvi.model.SCANVI.load_query_data(q,str(modeldir),accelerator='cpu')
 model.train(max_epochs=500,accelerator='cpu',devices=1,early_stopping=True,early_stopping_monitor='elbo_train',early_stopping_patience=10,early_stopping_min_delta=.001,plan_kwargs={'weight_decay':0.0},callbacks=[Progress()],enable_progress_bar=False)
 z=model.get_latent_representation();assert z.shape==(8174,30) and np.isfinite(z).all()
 model.save(str(O/'objects/HLCA_query_model'),overwrite=False)
 for key,v in model.history.items():
  if hasattr(v,'to_csv'):v.to_csv(O/'tables'/('HLCA_training_'+key+'.tsv'),sep='\t')
 np.save(O/'objects/HLCA_query_latent.npy',z)
 status('READING_REFERENCE_CORE')
 labels=[f'ann_level_{i}' for i in range(1,6)]+['ann_finest_level'];refcols={}
 with h5py.File(A/'HLCA_full_v1.1_emb.h5ad','r') as h:
  g=h['obs/core_or_extension'];cats=g['categories'].asstr()[:];mask=g['codes'][:]==list(cats).index('core')
  rx=h['X'][:][mask].astype(np.float32)
  for label in labels+['lung_condition','dataset','donor_id']:
   g=h['obs'][label]
   if isinstance(g,h5py.Group):
    cats=g['categories'].asstr()[:];codes=g['codes'][:][mask];refcols[label]=np.array([cats[i] if i>=0 else 'Unknown' for i in codes],dtype=object)
   else:refcols[label]=g.asstr()[:][mask]
 assert rx.shape[1]==30 and np.isfinite(rx).all()
 ref=pd.DataFrame(refcols);ref.groupby(['ann_finest_level','lung_condition']).size().reset_index(name='cells').to_csv(O/'tables/HLCA_reference_core_composition.tsv',sep='\t',index=False)
 status('REFERENCE_KNN',reference_cells=len(rx),query_cells=len(z))
 knn=KNeighborsTransformer(n_neighbors=50,mode='distance',algorithm='brute',metric='euclidean',n_jobs=4).fit(rx)
 ds=[];ixs=[]
 for start in range(0,len(z),256):
  with threadpool_limits(limits=4):dist,ix=knn.kneighbors(z[start:start+256])
  ds.append(dist);ixs.append(ix);status('REFERENCE_KNN',queried=min(start+256,len(z)),total=len(z))
 distances=np.vstack(ds);indices=np.vstack(ixs);sd=np.std(distances,axis=1)
 scale=np.divide(2.,sd,out=np.full_like(sd,np.inf),where=sd!=0)**2
 affin=np.exp(-distances/scale[:,None]);weights=affin/affin.sum(1,keepdims=True)
 assert np.isfinite(weights).all() and np.allclose(weights.sum(1),1)
 out=ev[['sample','dataset','donor_record','site','clinical_group','class_id','in_supported','dual_reference']].rename(columns={'class_id':'expanded_class_id'}).copy();out['supported_class_id']=su.class_id.reindex(out.index)
 probs=None
 for label in labels:
  values=ref[label].to_numpy();cats=np.array(sorted(set(values)),dtype=object)
  votes=np.column_stack([(weights*(values[indices]==cat)).sum(1) for cat in cats])
  best=votes.argmax(1);out[label+'_pred']=cats[best];out[label+'_uncertainty']=np.maximum(1-votes[np.arange(len(out)),best],0)
  if label=='ann_finest_level':
   probs=pd.DataFrame(votes,index=out.index,columns=cats);probs.to_csv(O/'tables/HLCA_finest_label_weights.tsv.gz',sep='\t',index_label='cell_id')
   for target in ['Adventitial fibroblasts','Alveolar fibroblasts','Pericytes','Smooth muscle','Myofibroblasts']:out[target.replace(' ','_')+'_weight']=probs[target] if target in probs else 0.
 out['nearest_reference_distance']=distances[:,0];out['mean50_reference_distance']=distances.mean(1)
 out.to_csv(O/'tables/HLCA_cell_predictions.tsv.gz',sep='\t',index_label='cell_id')
 summary=[]
 for layer,lab in [('supported',su),('expanded',ev)]:
  v=out.loc[lab.index].copy();v['class_id']=lab.class_id
  den=v.groupby(['dataset','donor_record','site','class_id'],dropna=False).size()
  for (study,donor,site,cls),g in v.groupby(['dataset','donor_record','site','class_id'],dropna=False):
   for pred,n in g.ann_finest_level_pred.value_counts().items():summary.append(dict(layer=layer,dataset=study,donor_record=donor,site=site,class_id=cls,reference_phenotype=pred,cells=int(n),total_candidate_cells=len(g),unfiltered_label_fraction=n/len(g),mean_uncertainty=g.loc[g.ann_finest_level_pred.eq(pred),'ann_finest_level_uncertainty'].mean(),no_confidence_cutoff=True))
  v.groupby(['dataset','donor_record','site','class_id'],dropna=False)[[c for c in v if c.endswith('_weight')]].mean().to_csv(O/'tables'/f'{layer}_HLCA_patient_phenotype_weights.tsv',sep='\t')
 pd.DataFrame(summary).to_csv(O/'tables/HLCA_patient_site_class_proportions.tsv',sep='\t',index=False)
 embed=ad.AnnData(z,obs=out.copy());clean(embed);embed.write_h5ad(O/'objects/HLCA_query_embedding_and_predictions.h5ad',compression='gzip')
 qa={'at':now(),'query_cells':len(out),'supported_cells':len(su),'reference_core_cells':len(rx),'reference_dimensions':30,'all_ids_aligned':out.index.equals(ev.index),'finite_latent':bool(np.isfinite(z).all()),'weights_sum_to_one':True,'uncertainties_in_0_1':bool(out[[c for c in out if c.endswith('_uncertainty')]].stack().between(0,1).all()),'confidence_cutoff':None,'scvi_version':scvi.__version__,'torch_version':torch.__version__,'no_primary_labels_replaced':True,'model_sha256':hashlib.sha256((modeldir/'model.pt').read_bytes()).hexdigest(),'official_mapping_not_CAF_confirmation':True}
 assert qa['all_ids_aligned'] and qa['uncertainties_in_0_1']
 writej(O/'delivery/HLCA_QA.json',qa);(O/'hlca.exit_code').write_text('0');(O/'HLCA_COMPLETE').write_text(now());status('MAPPING_COMPLETE_REVIEW_REQUIRED',**qa)
except Exception as e:
 traceback.print_exc();(O/'hlca.exit_code').write_text('1');status('FAILED',error=str(e));raise
