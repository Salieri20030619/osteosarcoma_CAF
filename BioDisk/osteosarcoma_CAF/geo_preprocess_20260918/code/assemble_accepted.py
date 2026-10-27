from pathlib import Path
import os,json,shutil
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ['MPLBACKEND']='Agg'
import numpy as np,pandas as pd,anndata as ad,matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from PIL import Image,ImageDraw
threadpool_limits(4)
P=Path.cwd();S=pd.read_csv(P/'delivery/GEO_sample_acceptance_and_QC.tsv',sep='\t')
D=pd.read_csv(P/'doublet_review/summary.tsv',sep='\t');assert len(D)==8 and (D.fraction<.25).all()
assert (P/'logs/doublet.exit_code').read_text().strip()=='0' and (P/'doublet_review/COMPLETE').exists()
root=P/'accepted';root.mkdir(exist_ok=True)
DOT=['CD3D','NKG7','CD79A','MZB1','LYZ','LST1','ACP5','CTSK','PECAM1','VWF','EPCAM','SFTPC','RGS5','NOTCH3','DCN','LUM','PDGFRA','FAP','RUNX2','ALPL','IBSP','ACTA2','TAGLN']
summary=[];changed=[];validations=[];cluster_counts=[]
for row in S.to_dict('records'):
 sm=row['sample'];d=root/sm;d.mkdir(exist_ok=True);a=ad.read_h5ad(P/'results'/sm/'annotated.h5ad')
 ev=pd.read_csv(P/'review'/(sm+'_cell_evidence.tsv.gz'),sep='\t',index_col='cell_id',dtype={'cluster':str})
 assert ev.index.tolist()==a.obs_names.tolist()
 a.obs['lineage_reviewed_candidate']=ev.lineage_review_proposed.to_numpy()
 for col in ['ref_eligible_immune','ref_eligible_endo','cnv_target_expanded','mixed_lineage_flag']:a.obs[col]=ev[col].to_numpy()
 a.obs['doublet_review_exclude']=False;a.obs['scDblFinder_score']=np.nan;a.obs['scDblFinder_class']='not_run_Scrublet_gate_passed'
 method='Scrublet_valid_automatic_threshold'
 if sm in set(D['sample']):
  calls=pd.read_csv(P/'doublet_review'/sm/'calls.tsv',sep='\t',index_col='cell_id')
  assert calls.index.is_unique and calls.index.tolist()==a.obs_names.tolist()
  a.obs['scDblFinder_score']=calls.scDblFinder_score.to_numpy();a.obs['scDblFinder_class']=calls.scDblFinder_class.to_numpy()
  a.obs['doublet_review_exclude']=calls.scDblFinder_class.eq('doublet').to_numpy();method='scDblFinder_after_Scrublet_threshold_review'
  fig,axs=plt.subplots(1,2,figsize=(10,4))
  axs[0].hist(a.obs.scDblFinder_score,bins=50,color='#27699b');axs[0].set(xlabel='scDblFinder score',ylabel='Cells',title=sm.split('__')[1])
  axs[1].scatter(a.obs.doublet_score,a.obs.scDblFinder_score,c=a.obs.doublet_review_exclude.map({False:'#27699b',True:'#bd483b'}),s=3)
  axs[1].set(xlabel='Scrublet score',ylabel='scDblFinder score',title='Red: alternate-method doublet call')
  fig.tight_layout();fig.savefig(d/'Doublet_review.png',dpi=150);plt.close(fig)
 excluded=int(a.obs.doublet_review_exclude.sum());before=a.n_obs
 a.obs['doublet_final_method']=method;a.obs['analysis_ready_keep']=~a.obs.doublet_review_exclude
 a.obs.to_csv(d/'all_preprocessed_cells_review.tsv.gz',sep='\t',index_label='cell_id')
 b=a[a.obs.analysis_ready_keep].copy()
 assert b.n_obs==before-excluded and (b.layers['counts']!=a.layers['counts'][a.obs.analysis_ready_keep.to_numpy()]).nnz==0
 b.obsm['X_umap_initial']=b.obsm.pop('X_umap')
 if 'X_pca' in b.obsm:b.obsm['X_pca_initial']=b.obsm.pop('X_pca')
 for key in list(b.obsp.keys()):del b.obsp[key]
 b.uns.pop('neighbors',None)
 b.uns['annotation_status']='broad_marker_reviewed_candidates_CNV_pending'
 b.uns['embedding_status']='initial_QC_embedding_coordinates_subset_after_doublet_review; no post-exclusion reclustering'
 b.uns['CAF_status']='NOT_CONFIRMED'
 b.write_h5ad(d/'analysis_ready.h5ad',compression='gzip')
 ck=ad.read_h5ad(d/'analysis_ready.h5ad',backed='r');assert ck.n_obs==b.n_obs and ck.obs_names.equals(b.obs_names);ck.file.close()
 coords=b.obsm['X_umap_initial'];labs=b.obs.lineage_reviewed_candidate.astype(str)
 fig,ax=plt.subplots(figsize=(8,6));palette=plt.get_cmap('tab20');classes=sorted(labs.unique())
 for i,label in enumerate(classes):
  mask=(labs==label).to_numpy();ax.scatter(coords[mask,0],coords[mask,1],s=3,color=palette(i%20),label=label,rasterized=True)
 ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',markerscale=3,frameon=False);ax.set(xlabel='Initial UMAP1',ylabel='Initial UMAP2',title=sm.split('__')[1]+' | marker-reviewed candidates')
 fig.tight_layout();fig.savefig(d/'UMAP.png',dpi=160,bbox_inches='tight');plt.close(fig)
 logx=b.X.tocsc();cnt=b.layers['counts'].tocsc();sym=b.var.symbol.astype(str).to_numpy()
 vals=[];detect=[]
 for gene in DOT:
  ix=np.flatnonzero(sym==gene)
  vals.append(np.asarray(logx[:,ix].sum(1)).ravel());detect.append(np.asarray(cnt[:,ix].sum(1)).ravel()>0)
 vals=np.array(vals).T;detect=np.array(detect).T
 clusters=sorted(b.obs.cluster.astype(str).unique(),key=int);means=[];pcts=[];labels=[]
 for cl in clusters:
  mask=(b.obs.cluster.astype(str)==cl).to_numpy();lab=labs[mask].unique();assert len(lab)==1
  means.append(vals[mask].mean(0));pcts.append(detect[mask].mean(0)*100);labels.append('C'+cl+' '+lab[0])
  cluster_counts.append(dict(sample=sm,cluster=cl,cells_after_doublet_review=int(mask.sum()),lineage_reviewed_candidate=lab[0]))
 means=np.array(means);pcts=np.array(pcts);sd=means.std(0);sd[sd==0]=1;z=(means-means.mean(0))/sd
 fig,ax=plt.subplots(figsize=(12,max(3,len(clusters)*.28+1.5)));xx,yy=np.meshgrid(np.arange(len(DOT)),np.arange(len(clusters)))
 pts=ax.scatter(xx.ravel(),yy.ravel(),s=pcts.ravel()*1.2,c=z.ravel(),cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
 ax.set_xticks(range(len(DOT)),DOT,rotation=90);ax.set_yticks(range(len(clusters)),labels);ax.invert_yaxis();ax.set_title(sm.split('__')[1]+' | reviewed candidates, CNV pending');ax.set_xlabel('Gene');ax.set_ylabel('Initial cluster')
 fig.colorbar(pts,ax=ax,label='Mean z');fig.tight_layout();fig.savefig(d/'Markers.png',dpi=160,bbox_inches='tight');plt.close(fig)
 shutil.copy2(P/'results'/sm/'QC.png',d/'QC_initial.png')
 summary.append(dict(sample=sm,analysis_ready_cells=b.n_obs,extra_doublets_excluded=excluded,final_doublet_method=method,
  fibroblast_like_after_doublet_review=int((labs=='FibLike').sum()),immune_reference_after_doublet_review=int(b.obs.ref_eligible_immune.sum()),endothelial_reference_after_doublet_review=int(b.obs.ref_eligible_endo.sum())))
 validations.append(dict(sample=sm,subset_count_conservation=True,saved_h5ad_cell_ids_and_shape_pass=True,before=before,excluded=excluded,after=b.n_obs))
 print('ACCEPTED_OBJECT',sm,b.n_obs,excluded,flush=True)
 del a,b,ck,logx,cnt
out=S.merge(pd.DataFrame(summary),on='sample',validate='one_to_one');out.to_csv(P/'delivery/GEO_sample_acceptance_and_QC.tsv',sep='\t',index=False)
pd.DataFrame(cluster_counts).to_csv(P/'delivery/GEO_cluster_counts_after_doublet_review.tsv',sep='\t',index=False)
(P/'review/accepted_output_integrity.json').write_text(json.dumps(validations,indent=2))
for kind in ['UMAP','Markers']:
 for i in range(0,len(out),4):
  batch=out.iloc[i:i+4];canvas=Image.new('RGB',(2200,1000*len(batch)),'white');draw=ImageDraw.Draw(canvas)
  for j,sm in enumerate(batch['sample']):
   im=Image.open(root/sm/(kind+'.png')).convert('RGB');im.thumbnail((2150,950));canvas.paste(im,(10,j*1000+25));draw.text((10,j*1000+5),sm,fill='black')
  canvas.save(P/'review'/('QA_accepted_'+kind+'_'+str(i//4+1)+'.png'))
(P/'accepted/COMPLETE').write_text('GEO_QC_AND_DOUBLETS_REVIEW_ASSEMBLED\n')
print(out[['sample','analysis_ready_cells','extra_doublets_excluded','fibroblast_like_after_doublet_review']].to_json(orient='records',indent=2),flush=True)
