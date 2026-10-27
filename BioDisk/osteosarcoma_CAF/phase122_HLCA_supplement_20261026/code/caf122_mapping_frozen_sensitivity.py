from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ['TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD']='1'
import numpy as np,pandas as pd,anndata as ad,json,h5py,scvi,torch
from lightning.pytorch import Callback
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';scvi.settings.seed=20261026;torch.set_num_threads(4)
q=ad.read_h5ad(O/'objects/query_model_counts.h5ad');trace=[]
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  rec=dict(additional_epoch=int(trainer.current_epoch)+1,metrics={k:float(v) for k,v in trainer.callback_metrics.items() if v.numel()==1});trace.append(rec)
  if len(trace)%10==0:print('TRAIN',len(trace),rec['metrics'],flush=True)
  (O/'resources/MAPPING_FROZEN_EXTENDED_TRAINING.json').write_text(json.dumps(trace),encoding='utf-8')
model=scvi.model.SCANVI.load_query_data(q,str(R/'phase32_lung_sensitivity_20260926/resources/HLCA_converted'),accelerator='cpu');model.module.load_state_dict({k:v for k,v in torch.load(O/'objects/query_model/model.pt',map_location='cpu',weights_only=False)['model_state_dict'].items() if k!='pyro_param_store'});frozen={k:v.detach().clone() for k,v in model.module.named_parameters() if not v.requires_grad};protected=model.module.z_encoder.encoder.fc_layers[0][0].weight[:,:2000].detach().clone();model.train(max_epochs=300,accelerator='cpu',devices=1,early_stopping=True,early_stopping_monitor='elbo_train',early_stopping_patience=10,early_stopping_min_delta=.001,plan_kwargs={'weight_decay':0.0},callbacks=[Progress()],enable_progress_bar=False,logger=False);freeze_max=max(float((dict(model.module.named_parameters())[k]-v).abs().max()) for k,v in frozen.items());protected_max=float((model.module.z_encoder.encoder.fc_layers[0][0].weight[:,:2000]-protected).abs().max());assert freeze_max==0 and protected_max==0;(O/'resources/FROZEN_RESUME_PARAMETER_CHECK.json').write_text(json.dumps(dict(status='PASS',frozen_parameters=len(frozen),frozen_max_change=freeze_max,encoder_gene_weight_max_change=protected_max)),encoding='utf-8');z=model.get_latent_representation();assert np.isfinite(z).all();np.save(O/'objects/query_latent_frozen_extended.npy',z);model.save(str(O/'objects/query_model_frozen_extended'),overwrite=False)
with h5py.File(R/'phase25_complete_nonPDX_assets_20260921/references/HLCA/HLCA_full_v1.1_emb.h5ad') as h:
 def col(k):
  g=h['obs'][k];return g['categories'].asstr()[:][g['codes'][:]] if isinstance(g,h5py.Group) else g.asstr()[:]
 core=col('core_or_extension')=='core';fine=col('ann_finest_level')[core];rx=h['X'][:][core].astype(np.float32)
def broad(v):
 if 'fibroblast' in v.lower():return 'Fibroblast'
 if v=='Pericytes':return 'Pericyte'
 if v.startswith(('Smooth muscle','SM activated')):return 'Smooth_muscle'
 return 'Other_lineage'
bc=np.array([broad(v) for v in fine]);cats=sorted(set(bc));knn=NearestNeighbors(n_neighbors=50,algorithm='brute',metric='euclidean',n_jobs=4).fit(rx);pred=[];finepred=[];allvotes=[]
for start in range(0,len(z),128):
 with threadpool_limits(limits=4):dist,ix=knn.kneighbors(z[start:start+128])
 sd=dist.std(1);scale=np.divide(2.,sd,out=np.full_like(sd,np.inf),where=sd!=0)**2;w=np.exp(-dist/scale[:,None]);w/=w.sum(1,keepdims=True);votes=np.column_stack([(w*(bc[ix]==c)).sum(1) for c in cats]);pred.extend(np.array(cats)[votes.argmax(1)]);allvotes.extend(votes)
 old=pd.read_csv(O/'tables/HLCA_current_cell_predictions.tsv.gz',sep='\t') if start==0 else old
new=old.copy();new['baseline_broad_pred']=old.broad_pred;new['broad_pred']=pred
for j,c in enumerate(cats):new['extended_vote_'+c]=np.array(allvotes)[:,j]
new.to_csv(O/'tables/HLCA_frozen_extended_cell_predictions.tsv.gz',sep='\t',index=False);summ=[]
for node in ['ALL','Lung_04','Bone_07','Bone_07_s1','Lung_04_s2']:
 g=new if node=='ALL' else new[new['leaf' if '_s' in node else 'parent'].eq(node)]
 summ.append(dict(node_id=node,cells=len(g),broad_agreement=float(g.broad_pred.eq(g.baseline_broad_pred).mean()),baseline_types=g.baseline_broad_pred.value_counts().to_dict(),extended_types=g.broad_pred.value_counts().to_dict()))
(O/'resources/MAPPING_FROZEN_EXTENDED_REVIEW.json').write_text(json.dumps(dict(status='PASS_COMPUTATION',baseline_epochs=200,additional_epochs=len(trace),results=summ),indent=2),encoding='utf-8');print(json.dumps(summ),flush=True)
