from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ['TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD']='1'
import numpy as np,pandas as pd,anndata as ad,json,datetime,hashlib,h5py,gc,scvi,torch
from sklearn.neighbors import NearestNeighbors
from lightning.pytorch import Callback
from threadpoolctl import threadpool_limits
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';old=R/'phase32_lung_sensitivity_20260926';refpath=R/'phase25_complete_nonPDX_assets_20260921/references/HLCA/HLCA_full_v1.1_emb.h5ad';modeldir=old/'resources/HLCA_converted'
def dump(n,x):(O/n).write_text(json.dumps(x,indent=2,default=str),encoding='utf-8')
def event(stage,**kw):
 r=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**kw);dump('MAPPING_STATUS.json',r)
 with (O/'logs/mapping_events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(r)+'\n')
 print(json.dumps(r),flush=True)
assert (O/'resources/FEATURE_ID_RECONCILIATION.json').exists()
scvi.settings.seed=20261026;torch.set_num_threads(4);q=ad.read_h5ad(O/'objects/query_model_counts.h5ad');assert q.n_obs==6423 and q.n_vars==2000
dump('resources/MAPPING_PROTOCOL.json',dict(method='HLCA official pretrained scANVI/scArches query adaptation; official50NN affinity voting',reference_model=str(modeldir),model_sha256=hashlib.sha256((modeldir/'model.pt').read_bytes()).hexdigest(),seed=20261026,query_batch='source dataset; no donor batch correction',epochs=200,early_stopping='elbo_train;patience10;min_delta0.001',reference_core_only_for_identity=True,missing_gene_policy='Model padding only; all feature coverage preserved, no missing-as-observed-zero inference',validation='Core pseudoquery with same-study donor neighbours removed; descriptive mapping sanity check, not independent disease validation'))
if not (O/'objects/query_latent.npy').exists():
 if not (O/'resources/MAPPING_PILOT_PASS.json').exists():
  ix=np.concatenate([g.sample(n=min(100,len(g)),random_state=20261026).index for _,g in q.obs.groupby('source',observed=True)]);pilot=q[ix].copy();scvi.model.SCANVI.prepare_query_anndata(pilot,str(modeldir));p=scvi.model.SCANVI.load_query_data(pilot,str(modeldir),accelerator='cpu');p.train(max_epochs=2,accelerator='cpu',devices=1,plan_kwargs={'weight_decay':0.0},enable_progress_bar=False,logger=False);zz=p.get_latent_representation();assert zz.shape==(len(ix),30) and np.isfinite(zz).all();dump('resources/MAPPING_PILOT_PASS.json',dict(status='PASS',cells=len(ix),finite_latent=True));del p,pilot;gc.collect()
 class Progress(Callback):
  def on_train_epoch_end(self,trainer,pl_module):
   if (trainer.current_epoch+1)%10==0:event('training',epoch=int(trainer.current_epoch)+1)
 event('training_start',cells=q.n_obs);scvi.model.SCANVI.prepare_query_anndata(q,str(modeldir));model=scvi.model.SCANVI.load_query_data(q,str(modeldir),accelerator='cpu');model.train(max_epochs=200,accelerator='cpu',devices=1,early_stopping=True,early_stopping_monitor='elbo_train',early_stopping_patience=10,early_stopping_min_delta=.001,plan_kwargs={'weight_decay':0.0},callbacks=[Progress()],enable_progress_bar=False,logger=False);z=model.get_latent_representation();assert z.shape==(6423,30) and np.isfinite(z).all();model.save(str(O/'objects/query_model'),overwrite=False);np.save(O/'objects/query_latent.npy',z)
 for k,v in model.history.items():v.to_csv(O/'tables'/('mapping_history_'+k+'.tsv'),sep='\t')
 del model;gc.collect()
else:z=np.load(O/'objects/query_latent.npy')
event('reading_reference')
with h5py.File(refpath) as h:
 def values(k):
  g=h['obs'][k]
  if isinstance(g,h5py.Group):
   cats=g['categories'].asstr()[:];codes=g['codes'][:];v=np.asarray(cats,dtype=object)[np.maximum(codes,0)];v[codes<0]='Unknown';return v
  return g.asstr()[:]
 columns=['core_or_extension','ann_finest_level','ann_level_1','ann_level_2','ann_level_3','ann_level_4','ann_level_5','lung_condition','dataset','study','donor_id','sample'];ref=pd.DataFrame({k:values(k) for k in columns});ref['cell_id']=values('_index');latent=h['X'][:].astype(np.float32)
core=ref.core_or_extension.eq('core').to_numpy();rx=latent[core];rm=ref[core].reset_index(drop=True);rm['donor_key']=rm.study.astype(str)+'::'+rm.donor_id.astype(str)
stroma=['Adventitial fibroblasts','Alveolar fibroblasts','Myofibroblasts','Peribronchial fibroblasts','Subpleural fibroblasts','Pericytes','Smooth muscle','Smooth muscle FAM83D+','SM activated stress response'];keep=ref.ann_finest_level.isin(stroma)&ref.lung_condition.isin(['Healthy','Healthy (tumor adjacent)','IPF']);ref[keep].to_csv(O/'tables/reference_stromal_metadata.tsv.gz',sep='\t',index=False);np.savez_compressed(O/'objects/reference_stromal_latent.npz',X=latent[keep]);rm.groupby(['ann_finest_level','study','donor_id']).size().rename('cells').reset_index().to_csv(O/'tables/reference_core_composition.tsv',sep='\t',index=False);del latent;gc.collect()
knn=NearestNeighbors(n_neighbors=50,algorithm='brute',metric='euclidean',n_jobs=4).fit(rx);distances=[];indices=[]
for start in range(0,len(z),128):
 with threadpool_limits(limits=4):dist,ix=knn.kneighbors(z[start:start+128])
 distances.append(dist);indices.append(ix)
 if start%1024==0:event('knn',cells=min(start+128,len(z)),total=len(z))
dist=np.vstack(distances);ix=np.vstack(indices);sd=dist.std(1);scale=np.divide(2.,sd,out=np.full_like(sd,np.inf),where=sd!=0)**2;w=np.exp(-dist/scale[:,None]);w/=w.sum(1,keepdims=True);out=q.obs[['cell_id','source','donor','site','parent','leaf','fibro_supported']].reset_index(drop=True).copy()
def broad(v):
 if 'fibroblast' in v.lower():return 'Fibroblast'
 if v=='Pericytes':return 'Pericyte'
 if v.startswith(('Smooth muscle','SM activated')):return 'Smooth_muscle'
 return 'Other_lineage'
fine=rm.ann_finest_level.to_numpy();bc=np.array([broad(x) for x in fine]);allvotes={}
for field,vals in [('finest',fine),('broad',bc)]:
 cats=sorted(set(vals));votes=np.column_stack([(w*(vals[ix]==a)).sum(1) for a in cats]);best=votes.argmax(1);out[field+'_pred']=np.array(cats)[best];out[field+'_uncertainty']=1-votes[np.arange(len(out)),best]
 for j,c in enumerate(cats):allvotes[field+'__'+c]=votes[:,j]
out['nearest_distance']=dist[:,0];out['mean50_distance']=dist.mean(1);assert np.isfinite(out.nearest_distance).all();out.to_csv(O/'tables/HLCA_current_cell_predictions.tsv.gz',sep='\t',index=False);votes=pd.concat([out[['cell_id']],pd.DataFrame(allvotes)],axis=1);votes.to_csv(O/'tables/HLCA_current_label_weights.tsv.gz',sep='\t',index=False)
nodes=['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2'];summary=[]
for node in nodes:
 col='leaf' if '_s' in node else 'parent';g=out[out[col].eq(node)]
 for (src,donor),a in g.groupby(['source','donor']):
  for typ in sorted(set(bc)):summary.append(dict(node_id=node,source=src,donor=donor,cells=len(a),identity=typ,predicted_fraction=float(a.broad_pred.eq(typ).mean()),mean_vote=float(votes.loc[a.index,'broad__'+typ].mean()),mean_uncertainty=float(a.broad_uncertainty.mean()),mean_distance=float(a.mean50_distance.mean())))
pd.DataFrame(summary).to_csv(O/'tables/candidate_identity_patient_summary.tsv',sep='\t',index=False)
# Neighbour exclusion validation keeps reference embeddings fixed and removes same donor from vote.
vi=rm.groupby('ann_finest_level',group_keys=False).sample(frac=.005,random_state=20261026).index.to_numpy();vi=vi[:3000];validation=[];knn.set_params(n_neighbors=200)
for start in range(0,len(vi),128):
 jj=vi[start:start+128]
 with threadpool_limits(limits=4):dd,ii=knn.kneighbors(rx[jj])
 for row,j in enumerate(jj):
  good=rm.donor_key.to_numpy()[ii[row]]!=rm.donor_key.iloc[j];sel=ii[row][good][:50]
  if not len(sel):validation.append(dict(index=int(j),label=fine[j],pred='Unassessable',available=0));continue
  pred=pd.Series(bc[sel]).value_counts().index[0];validation.append(dict(index=int(j),label=bc[j],pred=pred,available=len(sel),distance=float(dd[row][good][0])))
val=pd.DataFrame(validation);val.to_csv(O/'tables/reference_donor_exclusion_validation.tsv',sep='\t',index=False)
oldpred=pd.read_csv(old/'tables/HLCA_cell_predictions.tsv.gz',sep='\t');overlap=out.merge(oldpred[['cell_id','ann_finest_level_pred']],on='cell_id');overlap.to_csv(O/'tables/historical_mapping_overlap.tsv.gz',sep='\t',index=False)
dump('MAPPING_COMPLETE.json',dict(status='PASS_COMPUTATION',query_cells=len(out),core_reference_cells=len(rx),candidate_counts={n:int(out['leaf' if '_s' in n else 'parent'].eq(n).sum()) for n in nodes},reference_validation_cells=len(val),reference_validation_accuracy=float(val[val.available.gt(0)].eval('label == pred').mean()),historical_overlap_cells=len(overlap),historical_finest_agreement=float((overlap.finest_pred==overlap.ann_finest_level_pred).mean()),interpretation='reference similarity only; uncertainties and gene coverage retained'))
event('mapping_complete')
