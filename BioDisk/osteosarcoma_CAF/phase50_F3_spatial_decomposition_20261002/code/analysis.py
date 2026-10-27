from pathlib import Path
import sys,os,json,hashlib,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase46_spatial_functional_20261001';O=P/'phase50_F3_spatial_decomposition_20261002'
sys.path.insert(0,str(S/'code'))
from common import feature_map,corr,test_table,bh
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pr=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text())
z=pd.read_csv(S/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t',float_precision='round_trip',keep_default_na=False,na_values=[''],low_memory=False)
rows=[];checks=[];allspot=[];mapping=[]
for sm,t in z.groupby('sample',sort=True):
 t=t.copy();f=P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad'
 a=ad.read_h5ad(f);assert a.obs_names.astype(str).tolist()==t.spot_id.astype(str).tolist()
 gm=feature_map(a);idx=[gm[g] for g in pr['genes']]
 x=a.X.tocsr();lib=np.asarray(x.sum(axis=1)).ravel();assert np.array_equal(lib,t.library_UMI.to_numpy())
 raw=x[:,idx].toarray().astype(np.float64);log=np.log1p(raw/np.where(lib>0,lib,np.nan)[:,None]*1e4)
 for j,g in enumerate(pr['genes']):
  t[g]=log[:,j];mapping.append({'sample':sm,'gene':g,'feature_id':str(a.var_names[idx[j]]),'feature_index':idx[j],'nonzero_spots':int((raw[:,j]>0).sum())})
 ref=log.mean(axis=1);score_error=float(np.nanmax(np.abs(ref-t.F3_SELECTED_REMODELING.to_numpy(float))));assert np.allclose(ref,t.F3_SELECTED_REMODELING.to_numpy(float),atol=1e-12,rtol=0,equal_nan=True)
 adj=np.load(S/'objects'/f'{sm}_adjacency.npz');graph=sparse.csr_matrix((np.ones(len(adj['row'])),(adj['row'],adj['col'])),shape=tuple(adj['shape']));degree=np.asarray(graph.sum(1)).ravel()
 ng=np.divide(graph@t.CD8_RNA.to_numpy(float),degree,out=np.full(len(t),np.nan),where=degree>0)
 neighbor_error=float(np.nanmax(np.abs(ng-t.neighbor_CD8_RNA.to_numpy(float))));assert neighbor_error<1e-12
 checks.append({'sample':sm,'spots':len(t),'full_score_max_error':score_error,'neighbor_score_max_error':neighbor_error,'raw_input_SHA256':hashlib.sha256(f.read_bytes()).hexdigest()})
 for name,genes in pr['signatures'].items():t[name]=t[genes].mean(axis=1)
 for model,cols in pr['models'].items():
  c=t[cols].to_numpy(float) if cols else None
  if cols:
   for j,col in enumerate(cols):
    if col=='library_UMI':c[:,j]=np.log1p(c[:,j])
  for name in pr['signatures']:
   for ep in pr['endpoints']:
    rows.append({'sample':sm,'patient_id':t.patient_id.iloc[0],'site_group':t.site_group.iloc[0],'event_month':t.event_month.iloc[0],'signature':name,'model':model,'endpoint':ep,'effect':corr(t[name],t[ep],c),'spots':len(t)})
 keep=['spot_id','sample','patient_id','site_group','event_month','library_UMI','x','y','cancer','fibroblast','CD8_RNA','neighbor_CD8_RNA']+list(pr['signatures'])
 allspot.append(t[keep]);print('SAMPLE_COMPLETE',sm,len(t),flush=True)
 del a,x
f=pd.DataFrame(rows);f.to_csv(O/'tables/sample_effects.tsv',sep='\t',index=False)
keys=['signature','model','endpoint']
ev=f.groupby(['patient_id','site_group','event_month']+keys,dropna=False).effect.mean().reset_index()
sites=ev.groupby(['patient_id','site_group']+keys,dropna=False).effect.mean().reset_index()
patients=sites.groupby(['patient_id']+keys,dropna=False).effect.mean().reset_index()
assert patients.patient_id.nunique()==19 and len(f)==26*7*5*2
patients.to_csv(O/'tables/patient_effects.tsv',sep='\t',index=False);sites.to_csv(O/'tables/patient_site_effects.tsv',sep='\t',index=False)
stats=test_table(patients,keys,'effect',lambda r:'phase50_decomposition_'+r.model)
stats['FDR']=np.nan
for model,ids in stats[stats.signature.ne('full_reference')].groupby('model').groups.items():
 assert len(ids)==12
 stats.loc[ids,'FDR']=bh(stats.loc[ids,'P'])
stats['analysis_role']=np.where(stats.signature.eq('full_reference'),'previous_signature_reference','post_primary_exploratory')
stats['original_reference_FDR']=np.nan;old=pd.read_csv(S/'tables/spatial_patient_statistics.tsv',sep='\t',float_precision='round_trip');sens=pd.read_csv(S/'tables/spatial_model_sensitivity_statistics.tsv',sep='\t',float_precision='round_trip')
refchecks=[]
for i,r in stats[stats.signature.eq('full_reference')].iterrows():
 prev=None
 if r.model=='full_adjusted':prev=old[(old.program=='F3_SELECTED_REMODELING')&(old.endpoint=={'CD8_RNA':'CD8_adjusted','neighbor_CD8_RNA':'neighbor_CD8_adjusted'}[r.endpoint])]
 elif r.model=='raw' and r.endpoint=='CD8_RNA':prev=old[(old.program=='F3_SELECTED_REMODELING')&(old.endpoint=='CD8_raw')]
 elif r.endpoint=='CD8_RNA':prev=sens[(sens.program=='F3_SELECTED_REMODELING')&(sens.model==r.model)]
 if prev is not None and len(prev)==1:
  q=prev.iloc[0];err=abs(r.mean_effect-q.mean_effect);perr=abs(r.P-q.P);assert err<1e-12 and perr<1e-12
  stats.loc[i,'original_reference_FDR']=q.FDR
  refchecks.append({'model':r.model,'endpoint':r.endpoint,'mean_error':err,'P_error':perr,'original_FDR':q.FDR})
stats.to_csv(O/'tables/statistics.tsv',sep='\t',index=False)
loo=[]
for key,t in patients[patients.model.eq('full_adjusted')].groupby(keys):
 for pid in t.patient_id:
  d=t.loc[t.patient_id.ne(pid),'effect'].to_numpy(float);d=d[np.isfinite(d)]
  loo.append(dict(zip(keys,key),omitted_patient=pid,N=len(d),mean_effect=np.mean(d),negative_patients=int((d<0).sum()),positive_patients=int((d>0).sum())))
pd.DataFrame(loo).to_csv(O/'tables/leave_one_patient_out.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'tables/input_and_reconstruction_QA.tsv',sep='\t',index=False)
pd.DataFrame(mapping).to_csv(O/'tables/gene_feature_mapping.tsv',sep='\t',index=False)
pd.DataFrame(refchecks).to_csv(O/'tables/original_reference_reproduction_QA.tsv',sep='\t',index=False)
pd.concat(allspot,ignore_index=True).to_csv(O/'tables/spot_scores_for_QA.tsv.gz',sep='\t',index=False)
names=list(pr['signatures']);labels=['COL1A1','LOX','MMP14','LOX + MMP14','COL1A1 + MMP14','COL1A1 + LOX','Full 3-gene reference']
fig,axs=plt.subplots(1,2,figsize=(12,6),sharey=True)
for ax,ep in zip(axs,pr['endpoints']):
 for j,name in enumerate(names):
  v=patients[(patients.signature==name)&(patients.model=='full_adjusted')&(patients.endpoint==ep)].sort_values('patient_id').effect
  ax.scatter(v,np.full(len(v),j),s=16,alpha=.55,color='#527a9b')
  r=stats[(stats.signature==name)&(stats.model=='full_adjusted')&(stats.endpoint==ep)].iloc[0]
  ax.scatter([r.mean_effect],[j],s=60,marker='D',color='#ae492c',zorder=4)
 ax.axvline(0,color='gray',lw=.8);ax.set_yticks(range(7),labels);ax.set_title('Local CD8 RNA' if ep=='CD8_RNA' else 'Adjacent CD8 RNA');ax.set_xlabel('Patient mean partial rank correlation');ax.grid(axis='x',alpha=.2)
axs[0].invert_yaxis();fig.suptitle('F3 signal decomposition | 19 patients | full adjustment');fig.tight_layout();fig.savefig(O/'figures/patient_decomposition.png',dpi=170);fig.savefig(O/'figures/patient_decomposition.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,ep in zip(axs,pr['endpoints']):
 w=stats[stats.endpoint.eq(ep)].pivot(index='signature',columns='model',values='mean_effect').reindex(index=names,columns=list(pr['models']))
 im=ax.imshow(w,cmap='RdBu_r',vmin=-.12,vmax=.12)
 for j,name in enumerate(names):
  for k,model in enumerate(pr['models']):
   r=stats[(stats.signature==name)&(stats.model==model)&(stats.endpoint==ep)].iloc[0]
   ax.text(k,j,f'{r.mean_effect:.3f}',ha='center',va='center',fontsize=9,color='black')
 ax.set_yticks(range(7),labels);ax.set_xticks(range(5),['Raw','Library + XY','Cancer + library + XY','Fibroblast + library + XY','Full'],rotation=40,ha='right');ax.set_title('Local CD8 RNA' if ep=='CD8_RNA' else 'Adjacent CD8 RNA')
fig.subplots_adjust(left=.17,right=.86,bottom=.25,top=.9,wspace=.14);cax=fig.add_axes([.90,.30,.018,.48]);fig.colorbar(im,cax=cax,label='Patient mean rank correlation');fig.savefig(O/'figures/model_dependence.png',dpi=170);fig.savefig(O/'figures/model_dependence.pdf');plt.close(fig)
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(checks),'patients':19,'spots':sum(x['spots'] for x in checks),'statistics':len(stats),'new_hypotheses':60,'reference_rows':10,'FDR_families':5,'hypotheses_per_family':12,'original_reference_tests_reproduced':len(refchecks),'all_original_mean_P_errors_below_1e_12':True,'input_raw_counts_used':True,'no_spot_P':True,'leave_one_out_is_descriptive':True,'protocol_SHA256':hashlib.sha256((O/'resources/PROTOCOL_FROZEN.json').read_bytes()).hexdigest()}
(O/'delivery/ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2));(O/'ANALYSIS_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
