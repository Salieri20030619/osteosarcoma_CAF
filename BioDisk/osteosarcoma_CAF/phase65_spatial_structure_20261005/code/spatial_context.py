from pathlib import Path
import os,sys,json,gzip,io,gc,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005';S=B/'phase46_spatial_functional_20261001'
sys.path.insert(0,str(S/'code'));from common import feature_map,corr,test_table,sha
cfg=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text());struct=cfg['structure'];contexts=cfg['contexts'];drops=cfg['drop_dominant']
definitions={'full':struct,'drop_dominant_pair':{k:[g for g in gs if g not in drops[k]] for k,gs in struct.items()}}
base=pd.read_csv(S/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t',usecols=['spot_id','sample','patient_id','site_group','event_month','library_UMI','x','y'])
assert len(base)==53309 and base.patient_id.nunique()==19
records=[];frames=[];coverage=[];inputs=[]
for sm,t in base.groupby('sample',sort=True):
 t=t.copy();path=B/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad'
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==t.spot_id.tolist()
 x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();assert np.array_equal(lib,t.library_UMI)
 t['valid_library']=lib>0
 gm=feature_map(a);allgenes=sorted(set().union(*map(set,struct.values()),*map(set,contexts.values())));present=[g for g in allgenes if g in gm]
 val=np.log1p(x[:,[gm[g] for g in present]].toarray()/np.where(lib>0,lib,np.nan)[:,None]*1e4);lookup={g:i for i,g in enumerate(present)}
 # Store gene-level inputs so all correlations/scores can be independently recalculated without reloading large objects.
 for g in present:t['gene__'+g]=val[:,lookup[g]]
 for group,defs in [('context',contexts)]+list(definitions.items()):
  for name,genes in defs.items():
   have=[g for g in genes if g in lookup];ok=len(have)>=3 and len(have)/len(genes)>=.8
   t[group+'__'+name]=val[:,[lookup[g] for g in have]].mean(1) if ok else np.nan
   coverage.append({'sample':sm,'group':group,'program':name,'expected':len(genes),'measured':len(have),'evaluable':ok,'missing':';'.join(set(genes)-set(have)),'nonzero_spots':int((t[group+'__'+name]>0).sum())})
 for version in definitions:
  for program in struct:
   for ctx in contexts:
    for model in ['raw','depth_xy','mutual_contexts']:
     control=None
     if model!='raw':
      cols=[np.log1p(lib),t.x,t.y]
      if model=='mutual_contexts':cols.extend([t['context__'+c] for c in contexts if c!=ctx])
      control=np.column_stack(cols)
     r=corr(t[version+'__'+program],t['context__'+ctx],control)
     records.append({'sample':sm,'patient_id':t.patient_id.iloc[0],'site_group':t.site_group.iloc[0],'event_month':t.event_month.iloc[0],'version':version,'program':program,'context':ctx,'model':model,'effect':r,'spots':len(t)})
 frames.append(t)
 # Display all specimens; no best-looking section selected as representative.
 imagefile=Path(a.uns['provenance']['lowres_image_path']);scale=float(a.uns['spatial_scalefactors']['tissue_lowres_scalef'])
 im=Image.open(io.BytesIO(gzip.open(imagefile,'rb').read())).convert('RGB');xy=t[['x','y']].to_numpy()*scale
 assert ((xy[:,0]>=0)&(xy[:,0]<im.width)&(xy[:,1]>=0)&(xy[:,1]<im.height)).all()
 fig,axes=plt.subplots(2,4,figsize=(15,8),layout='constrained')
 axes.flat[0].imshow(im);axes.flat[0].set_title('H&E')
 plotkeys=['full__'+k for k in struct]+['context__'+k for k in contexts]
 for ax,k in zip(list(axes.flat)[1:],plotkeys):
  ax.imshow(im);sc=ax.scatter(xy[:,0],xy[:,1],c=t[k],s=4,alpha=.8,cmap='viridis');ax.set_title(k.replace('full__','').replace('context__',''),fontsize=10);fig.colorbar(sc,ax=ax,shrink=.55,pad=.015)
 for ax in axes.flat:ax.set_axis_off();ax.set_xlim(0,im.width);ax.set_ylim(im.height,0)
 fig.suptitle(f'{sm} | {t.patient_id.iloc[0]} | {t.site_group.iloc[0]}\nRNA context panels; no independent tumor contour or physical collagen measurement.',fontsize=12)
 fig.savefig(O/'figures/maps'/f'{sm}.png',dpi=125);plt.close(fig)
 inputs.append({'sample':sm,'path':str(path),'sha256':sha(path),'image_path':str(imagefile),'image_sha256':sha(imagefile),'spots':len(t)})
 del a,x,val;gc.collect();print('SPATIAL',sm,len(t),flush=True)
effects=pd.DataFrame(records);effects.to_csv(O/'tables/spatial_sample_correlations.tsv',sep='\t',index=False)
pd.concat(frames,ignore_index=True).to_csv(O/'tables/spatial_gene_and_program_scores.tsv.gz',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'tables/spatial_coverage.tsv',sep='\t',index=False)
pd.DataFrame(inputs).to_csv(O/'resources/SPATIAL_INPUT_MANIFEST.tsv',sep='\t',index=False)
keys=['version','program','context','model'];allpat=[];allstats=[];allcontrast=[];contraststats=[]
for scope,sub in [('all_patients',effects),('primary_only',effects[effects.site_group=='primary_tumor'])]:
 events=sub.groupby(['patient_id','site_group','event_month']+keys,dropna=False).effect.mean().reset_index()
 sites=events.groupby(['patient_id','site_group']+keys,dropna=False).effect.mean().reset_index()
 pat=sites.groupby(['patient_id']+keys,dropna=False).effect.mean().reset_index();pat['cohort_scope']=scope;allpat.append(pat)
 stats=test_table(pat,keys,'effect',lambda r:scope+'|'+r.version+'|'+r.model+'|context12');stats['cohort_scope']=scope;allstats.append(stats)
 con=[]
 for k,d in pat.groupby(['patient_id','version','program','model']):
  v=d.set_index('context').effect
  for a,b in cfg['contrasts']:con.append(dict(zip(['patient_id','version','program','model'],k),contrast=a+'_minus_'+b,effect=v[a]-v[b]))
 con=pd.DataFrame(con);con['cohort_scope']=scope;allcontrast.append(con)
 cs=test_table(con,['version','program','model','contrast'],'effect',lambda r:scope+'|'+r.version+'|'+r.model+'|source_preference9');cs['cohort_scope']=scope;contraststats.append(cs)
patients=pd.concat(allpat,ignore_index=True);stats=pd.concat(allstats,ignore_index=True)
patients.to_csv(O/'tables/spatial_patient_correlations.tsv',sep='\t',index=False);stats.to_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t',index=False)
pd.concat(allcontrast,ignore_index=True).to_csv(O/'tables/source_preference_patient_effects.tsv',sep='\t',index=False)
pd.concat(contraststats,ignore_index=True).to_csv(O/'tables/source_preference_statistics.tsv',sep='\t',index=False)
rob=[]
for k,d in patients.groupby(['cohort_scope']+keys):
 v=d.effect.dropna().to_numpy();n=len(v);loo=(v.sum()-v)/(n-1) if n>1 else np.array([np.nan])
 rob.append(dict(zip(['cohort_scope']+keys,k),N=n,mean_effect=v.mean() if n else np.nan,loo_min=np.nanmin(loo) if n>1 else np.nan,loo_max=np.nanmax(loo) if n>1 else np.nan,all_loo_same_direction=bool(np.all(np.sign(loo)==np.sign(v.mean()))) if n>1 else False))
pd.DataFrame(rob).to_csv(O/'tables/spatial_patient_leave_one_out.tsv',sep='\t',index=False)
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,scope in zip(axes,['all_patients','primary_only']):
 sub=stats[(stats.cohort_scope==scope)&(stats.version=='full')&(stats.model=='depth_xy')]
 w=sub.pivot(index='program',columns='context',values='mean_effect')
 im=ax.imshow(w,vmin=-.6,vmax=.6,cmap='RdBu_r',aspect='auto')
 ax.set_xticks(range(len(w.columns)),[x.replace('_context','') for x in w.columns],rotation=30,ha='right');ax.set_yticks(range(len(w.index)),w.index);ax.set_title(scope)
 for i,p in enumerate(w.index):
  for j,c in enumerate(w.columns):
   r=sub[(sub.program==p)&(sub.context==c)].iloc[0];ax.text(j,i,f'{r.mean_effect:.2f}'+('*' if r.P<.05 else ''),ha='center',va='center',fontsize=10)
 fig.colorbar(im,ax=ax,shrink=.75)
fig.suptitle('Patient mean partial rank correlation (depth, x, y adjusted)\n* Patient P < 0.05; full family FDR in tables. RNA association only.',fontsize=12)
fig.savefig(O/'figures/source_context_summary.png',dpi=170);fig.savefig(O/'figures/source_context_summary.pdf');plt.close(fig)
(O/'SPATIAL_COMPLETE.json').write_text(json.dumps({'success':True,'slices':len(inputs),'spots':len(base),'patients':19,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print('SPATIAL_COMPLETE',len(stats),flush=True)
