from pathlib import Path
import json,itertools,hashlib
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase9_joint_20260919'
a=ad.read_h5ad(O/'objects/joint_CAF_integrated.h5ad');ob=a.obs.copy();ob['joint_cluster']=ob.joint_cluster.astype(str)
clusters=sorted(ob.joint_cluster.unique(),key=int);x=a.layers['counts'].tocsr();lib=ob.full_library_counts.to_numpy();sy=a.var.symbol.astype(str).to_numpy()
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
support=pd.read_csv(O/'tables/cluster_support.tsv',sep='\t');s=support[support.analysis=='joint_main'].copy().sort_values('cluster')
s['interpretation']=np.where(s.max_known_donor_fraction>.7,'donor_dominated',np.where(s.known_donor_records_ge20<3,'limited_donor_support','multi_donor_candidate'))
save(s,'joint_cluster_review.tsv')
save(ob.groupby(['joint_cluster','dataset','sample','patient_id','donor_known','site'],observed=True).size().rename('cells').reset_index(),'joint_cluster_sample_counts.tsv')
# Cell-pooled ranking is descriptive only; no cell-level inferential test.
rank=[]
for cl in clusters:
 m=ob.joint_cluster.to_numpy()==cl
 mean=np.asarray(a.X[m].mean(0)).ravel();rest=np.asarray(a.X[~m].mean(0)).ravel();det=np.asarray((x[m]>0).mean(0)).ravel();dr=np.asarray((x[~m]>0).mean(0)).ravel()
 d=pd.DataFrame({'cluster':cl,'ensembl_id':a.var_names,'symbol':sy,'mean_logCP10k_difference':mean-rest,'detection_cluster':det,'detection_other':dr})
 d['rank_score']=d.mean_logCP10k_difference*np.maximum(det-dr,0)
 d=d.sort_values('rank_score',ascending=False);d['descriptive_rank']=np.arange(1,len(d)+1);rank.append(d)
rank=pd.concat(rank);save(rank,'joint_descriptive_all_gene_markers.tsv.gz');save(rank.groupby('cluster',sort=False).head(30),'joint_descriptive_top30.tsv')
# Patient-paired pseudobulk within each study. Full raw library denominator.
results=[];pairs=[]
for ds in sorted(ob.dataset.unique()):
 if ds=='GSE250015':continue
 for cl in clusters:
  deltas=[]
  for donor in ob.loc[(ob.dataset==ds)&ob.donor_known,'donor_record'].unique():
   dm=(ob.donor_record.to_numpy()==donor);m=dm&(ob.joint_cluster.to_numpy()==cl);r=dm&~m
   if m.sum()<20 or r.sum()<20:continue
   v=np.log2(1+np.asarray(x[m].sum(0)).ravel()/lib[m].sum()*1e6)-np.log2(1+np.asarray(x[r].sum(0)).ravel()/lib[r].sum()*1e6)
   deltas.append(v);pairs.append({'dataset':ds,'cluster':cl,'donor_record':donor,'cluster_cells':int(m.sum()),'other_cells':int(r.sum())})
  if not deltas:continue
  z=np.vstack(deltas);n=len(z);p=np.full(a.n_vars,np.nan);q=p.copy()
  if n>=3:
   signs=np.array(list(itertools.product([-1,1],repeat=n)))
   null=np.abs(signs@z/n);stat=np.abs(z.mean(0));p=(null>=stat[None,:]-1e-12).mean(0);q=multipletests(p,method='fdr_bh')[1]
  results.append(pd.DataFrame({'dataset':ds,'cluster':cl,'ensembl_id':a.var_names,'symbol':sy,'paired_patients':n,'median_log2CPM_difference':np.median(z,axis=0),'positive_patients':(z>0).sum(0),'negative_patients':(z<0).sum(0),'P_exact_signflip':p,'FDR_BH':q,'evidence':'exploratory_exact_test' if n>=3 else 'descriptive_only'}))
res=pd.concat(results);save(res,'joint_within_study_patient_effects.tsv.gz');save(pd.DataFrame(pairs),'pseudobulk_patient_pairs.tsv')
# Match clusters by actual cell overlap, including subset sensitivity runs.
ov=[]
for path in sorted((O/'tables').glob('*_cell_labels.tsv.gz')):
 name=path.name.replace('_cell_labels.tsv.gz','')
 if name=='joint_main':continue
 b=pd.read_csv(path,sep='\t',index_col=0);main=ob.loc[b.index,'joint_cluster'];labels=b[name].astype(str)
 for cl in clusters:
  m=main==cl
  if not m.any():continue
  scores=[]
  for other in labels.unique():
   v=labels==other;inter=int((m&v).sum());scores.append((inter/int((m|v).sum()),other,inter,int(v.sum())))
  jac,other,inter,sz=max(scores)
  ov.append({'analysis':name,'joint_cluster':cl,'main_cells_in_subset':int(m.sum()),'best_cluster':other,'overlap_cells':inter,'best_cluster_cells':sz,'Jaccard':jac,'main_recall':inter/int(m.sum())})
save(pd.DataFrame(ov),'cluster_sensitivity_overlap.tsv')
plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
fig,axs=plt.subplots(2,3,figsize=(17,10))
for row,rep in enumerate(['X_umap_uncorrected','X_umap_joint_main']):
 for col,key in enumerate(['dataset','site','joint_cluster']):
  ax=axs[row,col];vals=ob[key].astype(str);cats=sorted(vals.unique(),key=int) if key=='joint_cluster' else sorted(vals.unique());palette=plt.get_cmap('tab20')
  for j,v in enumerate(cats):
   m=(vals==v).to_numpy();xy=a.obsm[rep][m];ax.scatter(xy[:,0],xy[:,1],s=3,c=[palette(j)],label=('J'+v if key=='joint_cluster' else v),rasterized=True,alpha=.8)
  ax.set_title(('Before correction' if row==0 else 'Harmony: study + technology')+' | '+key);ax.set_xticks([]);ax.set_yticks([]);ax.legend(markerscale=3,fontsize=8,loc='best',framealpha=.8)
fig.tight_layout();fig.savefig(O/'figures/01_joint_umap.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(16,6));
for ax,key in zip(axs,['sample_key','donor_record']):
 vals=ob[key].astype(str)
 for j,v in enumerate(sorted(vals.unique())):
  m=vals.eq(v).to_numpy();xy=a.obsm['X_umap_joint_main'][m];ax.scatter(xy[:,0],xy[:,1],s=3,color=plt.get_cmap('tab20')(j),label=v)
 ax.legend(fontsize=7,markerscale=3,loc='center left',bbox_to_anchor=(1,0.5));ax.set_title('Joint Harmony | '+key);ax.set_xticks([]);ax.set_yticks([])
fig.tight_layout();fig.savefig(O/'figures/02_sample_donor.png',dpi=180,bbox_inches='tight');plt.close(fig)
# Program heatmap uses standardized cluster mean, for description only.
progs=pd.read_csv(O/'tables/cluster_programs.tsv',sep='\t');t=progs[progs.analysis=='joint_main'].pivot(index='program',columns='cluster',values='mean_score').reindex(columns=[int(c) for c in clusters]);z=t.sub(t.mean(1),axis=0).div(t.std(1).replace(0,np.nan),axis=0).fillna(0)
fig,ax=plt.subplots(figsize=(13,10));im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(z)),z.index);ax.set_xticks(range(12),['J'+c for c in clusters]);ax.set_title('Program profiles | cluster means standardized across clusters');fig.colorbar(im,ax=ax,label='Z score of cluster mean');fig.tight_layout();fig.savefig(O/'figures/03_program_profiles.png',dpi=180);plt.close(fig)
ct=pd.crosstab(ob.joint_cluster,ob.dataset).reindex(clusters);ct=ct.div(ct.sum(1),axis=0)
t=pd.DataFrame(ov).pivot(index='analysis',columns='joint_cluster',values='Jaccard').reindex(columns=clusters)
fig,axs=plt.subplots(1,2,figsize=(16,6));ct.plot.bar(stacked=True,ax=axs[0],colormap='tab10');axs[0].set_xlabel('Joint cluster');axs[0].set_ylabel('Fraction of cells');axs[0].legend(fontsize=8);axs[0].set_title('Study composition of each cluster')
im=axs[1].imshow(t,cmap='viridis',vmin=0,vmax=1,aspect='auto');axs[1].set_yticks(range(len(t)),t.index);axs[1].set_xticks(range(12),['J'+c for c in clusters]);axs[1].set_title('Best cluster overlap in each sensitivity');fig.colorbar(im,ax=axs[1],label='Jaccard on shared cells');fig.tight_layout();fig.savefig(O/'figures/04_support_sensitivity.png',dpi=180);plt.close(fig)
# Separate raw marker panel, not Harmony-corrected gene values.
markers=['DCN','LUM','PDGFRA','FAP','THY1','ACTA2','TAGLN','MYL9','CXCL12','IL6','CCL2','HLA-DRA','CD74','COL1A1','COL3A1','POSTN','CTHRC1','MGP','ALPL','RUNX2','S100A4','PTPRC','PECAM1','EPCAM','RGS5','CSPG4']
markers=[g for g in markers if g in sy];means=[];dets=[]
for cl in clusters:
 m=ob.joint_cluster.to_numpy()==cl;v=[];d=[]
 for g in markers:
  ix=np.flatnonzero(sy==g);raw=np.asarray(x[m][:,ix].sum(1)).ravel();v.append(np.log1p(raw/lib[m]*1e4).mean());d.append((raw>0).mean())
 means.append(v);dets.append(d)
means=np.array(means);dets=np.array(dets);norm=(means-means.mean(0))/np.maximum(means.std(0),1e-9)
fig,ax=plt.subplots(figsize=(15,6));xx,yy=np.meshgrid(np.arange(len(markers)),np.arange(12));im=ax.scatter(xx.ravel(),yy.ravel(),s=dets.ravel()*180,c=norm.ravel(),cmap='RdBu_r',vmin=-2,vmax=2);ax.set_xticks(range(len(markers)),markers,rotation=60,ha='right');ax.set_yticks(range(12),['J'+c for c in clusters]);ax.invert_yaxis();ax.set_title('Raw RNA markers | dot area = detection fraction; color = standardized mean logCP10k');fig.colorbar(im,ax=ax,label='Z score');fig.tight_layout();fig.savefig(O/'figures/05_marker_panel.png',dpi=180);plt.close(fig)
# Verify original input hashes recorded prior to integration.
man=pd.read_csv(O/'resources/source_manifest.tsv',sep='\t');checks=[]
for row in man.to_dict('records'):
 path=Path(row.get('path',row.get('source_file','')));actual=hashlib.sha256(path.read_bytes()).hexdigest();assert actual==row['sha256'];checks.append(str(path))
assert len(ob)==4300 and int(s.cells.sum())==4300 and x.shape==(4300,15391)
review={'source_files_hash_verified':len(checks),'cell_count':4300,'cluster_count':12,'donor_dominated_clusters':s.loc[s.interpretation=='donor_dominated','cluster'].tolist(),'multi_donor_candidate_clusters':s.loc[s.interpretation=='multi_donor_candidate','cluster'].tolist(),'within_study_pseudobulk_comparisons':len(res[['dataset','cluster']].drop_duplicates()),'comparisons_with_ge3_paired_patients':len(res.loc[res.paired_patients>=3,['dataset','cluster']].drop_duplicates()),'minimum_paired_patients':int(res.paired_patients.min()),'maximum_paired_patients':int(res.paired_patients.max()),'FDR_lt_0_05_rows':int((res.FDR_BH<.05).sum()),'all_checks_passed':True}
(O/'REVIEW_QA.json').write_text(json.dumps(review,indent=2));print(json.dumps(review),flush=True);print(s.to_string(index=False));print(pd.DataFrame(pairs).groupby(['dataset','cluster']).size().to_string())
(O/'REVIEW_COMPLETE').write_text('REVIEW_CHECKS_PASSED\n')