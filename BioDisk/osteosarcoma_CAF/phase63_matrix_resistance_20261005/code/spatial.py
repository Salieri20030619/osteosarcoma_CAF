from pathlib import Path
import os,sys,json,gzip,io,gc
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005';S=P/'phase46_spatial_functional_20261001'
sys.path.insert(0,str(S/'code'))
from common import feature_map,corr,test_table,canon,sha
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
cfg=json.loads((O/'resources/PROGRAMS.json').read_text());sets=cfg['sets'];roles=cfg['roles'];assert (O/'FREEZE_COMPLETE').exists()
ann=pd.read_excel(P/'phase35_external_validation_prep_20260927/sources/GSE293065_supplement_dataset.xlsx',header=1);groups=json.loads((S/'resources/AUTHOR_COMPONENT_GROUPS.json').read_text());used=set(sum(groups.values(),[]));contrib=set()
for row in ann[ann.IC.isin(used)].itertuples():
 contrib.update(canon(g.strip()) for g in str(row.contrib_genes).split(',') if canon(g.strip()))
reduced={k:sorted(set(v)-contrib) for k,v in sets.items()}
(O/'resources/COMPONENT_GENE_OVERLAP.json').write_text(json.dumps({'component_contributor_genes':sorted(contrib),'removed_by_program':{k:sorted(set(v)&contrib) for k,v in sets.items()},'reduced_definitions':reduced},indent=2))
base=read(S/'tables/all_spot_RNA_and_components.tsv.gz',float_precision='round_trip');assert len(base)==53309 and base.patient_id.nunique()==19
out=[];spotframes=[];cov=[];inputs=[];subset=['collagen_biosynthesis','collagen_crosslinking','ECM_degradation'];figdir=O/'figures/spatial_maps';figdir.mkdir(parents=True,exist_ok=True)
endpoint_specs={
'tumor_local_raw':('cancer',[],'main_tumor_neighborhood'),
'tumor_local_adjusted':('cancer',['fibroblast','library_UMI','x','y'],'main_tumor_neighborhood'),
'tumor_neighbor_adjusted':('neighbor_cancer',['cancer','fibroblast','library_UMI','x','y'],'main_tumor_neighborhood'),
'fibroblast_local_adjusted':('fibroblast',['cancer','library_UMI','x','y'],'source_context'),
'CD8_local_raw':('CD8_RNA',[],'eggshell_immune_side'),
'CD8_local_adjusted':('CD8_RNA',['cancer','fibroblast','library_UMI','x','y'],'eggshell_immune_side'),
'CD8_neighbor_adjusted':('neighbor_CD8_RNA',['cancer','fibroblast','library_UMI','x','y'],'eggshell_immune_side')}
for sm,t in base.groupby('sample',sort=True):
 t=t.copy();path=P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad';a=ad.read_h5ad(path);assert a.obs_names.tolist()==t.spot_id.tolist();x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();assert np.array_equal(lib,t.library_UMI)
 gm=feature_map(a);allgenes=sorted(set().union(*map(set,sets.values())));present=[g for g in allgenes if g in gm];val=np.log1p(x[:,[gm[g] for g in present]].toarray()/np.where(lib>0,lib,np.nan)[:,None]*1e4);lookup={g:j for j,g in enumerate(present)}
 adj=np.load(S/'objects'/f'{sm}_adjacency.npz');graph=sparse.csr_matrix((np.ones(len(adj['row'])),(adj['row'],adj['col'])),shape=tuple(adj['shape']));deg=np.asarray(graph.sum(1)).ravel()
 for name in ['cancer','CD8_RNA']:
  vals=t[name].to_numpy(float);nv=np.divide(graph@vals,deg,out=np.full(len(t),np.nan),where=deg>0)
  if name=='CD8_RNA':assert np.allclose(nv,t.neighbor_CD8_RNA,equal_nan=True)
  else:t['neighbor_cancer']=nv
 controls={ep:(np.column_stack([np.log1p(t[c]) if c=='library_UMI' else t[c] for c in cs]) if cs else None) for ep,(target,cs,role) in endpoint_specs.items()}
 for version,definitions in [('full_program',sets),('minus_author_contributors',reduced)]:
  for name,gs in definitions.items():
   have=[g for g in gs if g in lookup];ok=len(have)>=min(3,len(sets[name])) and len(gs)>0 and len(have)/len(gs)>=.8;score=val[:,[lookup[g] for g in have]].mean(1) if ok else np.full(len(t),np.nan);t[version+'__'+name]=score
   cov.append({'sample':sm,'version':version,'program':name,'original_genes':len(sets[name]),'definition_genes':len(gs),'measured_genes':len(have),'evaluable':ok,'missing':';'.join(sorted(set(gs)-set(have)))})
   for ep,(target,cs,role) in endpoint_specs.items():
    effect=corr(score,t[target],controls[ep]) if ok else np.nan
    out.append({'sample':sm,'patient_id':t.patient_id.iloc[0],'site_group':t.site_group.iloc[0],'event_month':t.event_month.iloc[0],'version':version,'program':name,'endpoint':ep,'role':role,'effect':effect,'spots':len(t)})
 keep=['spot_id','sample','patient_id','site_group','event_month','library_UMI','x','y','cancer','fibroblast','neighbor_cancer','CD8_RNA','neighbor_CD8_RNA']+[v+'__'+n for v in ['full_program','minus_author_contributors'] for n in sets];spotframes.append(t[keep])
 imfile=Path(a.uns['provenance']['lowres_image_path']);scale=float(a.uns['spatial_scalefactors']['tissue_lowres_scalef']);im=Image.open(io.BytesIO(gzip.open(imfile,'rb').read())).convert('RGB');xy=t[['x','y']].to_numpy()*scale;assert ((xy[:,0]>=0)&(xy[:,0]<im.width)&(xy[:,1]>=0)&(xy[:,1]<im.height)).all()
 fig,axs=plt.subplots(2,3,figsize=(13,8));names=['cancer','fibroblast','CD8_RNA']+['full_program__'+n for n in subset];titles=['Tumor signed IC','Fibroblast signed IC','CD8 RNA','Collagen biosynthesis RNA','Collagen crosslinking RNA','ECM degradation RNA']
 for ax,name,title in zip(axs.ravel(),names,titles):
  ax.imshow(im);sc=ax.scatter(xy[:,0],xy[:,1],c=t[name],s=6,cmap='viridis',alpha=.75);ax.set_xlim(0,im.width);ax.set_ylim(im.height,0);ax.set_aspect('equal');ax.set_axis_off();ax.set_title(title);fig.colorbar(sc,ax=ax,shrink=.65,pad=.015)
 fig.suptitle(sm+' | '+str(t.patient_id.iloc[0])+' | '+str(t.site_group.iloc[0]),fontsize=14);fig.text(.015,.012,'RNA programs and signed components over tissue image. No independently annotated tumor boundary or physical collagen measurement.',fontsize=9);fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(figdir/f'{sm}.png',dpi=150);plt.close(fig)
 inputs.append({'sample':sm,'spots':len(t),'raw_path':str(path),'raw_sha256':sha(path),'image_path':str(imfile),'adjacency_edges':int(graph.nnz/2),'independent_tumor_mask':False})
 del a,x,val;gc.collect();print('SPATIAL_SAMPLE',sm,len(t),flush=True)
f=pd.DataFrame(out);f.to_csv(O/'tables/spatial_sample_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/spatial_program_coverage.tsv',sep='\t',index=False);pd.DataFrame(inputs).to_csv(O/'resources/spatial_input_manifest.tsv',sep='\t',index=False)
keys=['version','program','endpoint','role'];events=f.groupby(['patient_id','site_group','event_month']+keys,dropna=False).effect.mean().reset_index();sites=events.groupby(['patient_id','site_group']+keys,dropna=False).effect.mean().reset_index();patients=sites.groupby(['patient_id']+keys,dropna=False).effect.mean().reset_index()
patients.to_csv(O/'tables/spatial_patient_effects.tsv',sep='\t',index=False);sites.to_csv(O/'tables/spatial_patient_site_effects.tsv',sep='\t',index=False)
stats=test_table(patients,keys,'effect',lambda r:'spatial|'+r.version+'|'+r.endpoint+'|'+roles[r.program]);stats.to_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t',index=False)
site_stats=test_table(sites,['site_group']+keys,'effect',lambda r:'spatial_site|'+r.site_group+'|'+r.version+'|'+r.endpoint+'|'+roles[r.program]);site_stats.to_csv(O/'tables/spatial_site_statistics.tsv',sep='\t',index=False)
pd.concat(spotframes,ignore_index=True).to_csv(O/'tables/spatial_spot_scores.tsv.gz',sep='\t',index=False)
main=stats[(stats.version=='full_program')&stats.program.isin(list(sets)[:7])];wide=main.pivot(index='program',columns='endpoint',values='mean_effect');fig,ax=plt.subplots(figsize=(12,6));lim=max(.15,np.nanmax(np.abs(wide.to_numpy())));im=ax.imshow(wide,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set_yticks(range(len(wide)),wide.index.str.replace('_',' '));ax.set_xticks(range(len(wide.columns)),wide.columns.str.replace('_',' '),rotation=35,ha='right')
for i,p in enumerate(wide.index):
 for j,e in enumerate(wide.columns):
  rr=main[(main.program==p)&(main.endpoint==e)].iloc[0];mark='*' if rr.P<.05 else '';ax.text(j,i,f'{rr.mean_effect:.2f}{mark}',ha='center',va='center',fontsize=9)
fig.colorbar(im,ax=ax,label='Patient mean rank correlation');ax.set_title('Matrix mechanisms in spatial tissue | 19 patient records');fig.text(.02,.015,'* Patient P < 0.05, nominal only; complete FDR in tables. Signed IC values are not cell fractions. No resistance or shell geometry endpoint implied.',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(O/'figures/spatial_mechanism_summary.png',dpi=170);plt.close(fig)
(O/'SPATIAL_COMPLETE').write_text('SUCCESS\n');print('SPATIAL_COMPLETE',len(stats),len(inputs),flush=True)
