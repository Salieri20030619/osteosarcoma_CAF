from pathlib import Path
import json,hashlib,warnings
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';assert (O/'BATCH_MODELS_COMPLETE').exists()
a=ad.read_h5ad(O/'objects/unified_strict_integrated.h5ad');ob=a.obs.copy();ob['unified_cluster']=ob.unified_cluster.astype(str);clusters=sorted(ob.unified_cluster.unique(),key=int);x=a.layers['counts'];sy=a.var.symbol.astype(str).to_numpy()
def save(d,name):d.to_csv(O/'tables'/name,sep='\t',index=False)
changes=pd.read_csv(O/'tables/sample_gate_changes.tsv',sep='\t');save(changes.groupby('dataset')[['retained_cells','old_joint_cells','unified_strict','osteo_sensitive','competitor_sensitive','mural_reference','old_retained_in_new','old_removed','new_vs_old']].sum().reset_index(),'cohort_gate_changes.tsv')
support=pd.read_csv(O/'tables/cluster_support.tsv',sep='\t');s=support[support.analysis=='unified_primary'].copy();s['interpretation']=np.where(s.max_known_donor_fraction>.7,'donor_dominated',np.where(s.donor_records_ge20<3,'limited_support','multiple_records_candidate'));save(s,'primary_cluster_review.tsv')
# Raw marker rankings remain descriptive and are separate from batch-corrected coordinates.
ranks=[];programs=[]
for c,g in ob.groupby('unified_cluster',observed=True):
 m=ob.unified_cluster.eq(c).to_numpy();mean=np.asarray(a.X[m].mean(0)).ravel()-np.asarray(a.X[~m].mean(0)).ravel();det=np.asarray((x[m]>0).mean(0)).ravel();other=np.asarray((x[~m]>0).mean(0)).ravel();d=pd.DataFrame({'cluster':c,'ensembl_id':a.var_names,'symbol':sy,'mean_logCP10k_difference':mean,'detection_cluster':det,'detection_other':other,'rank_score':mean*np.maximum(det-other,0)}).sort_values('rank_score',ascending=False);d['descriptive_rank']=np.arange(1,len(d)+1);ranks.append(d)
 for name in [c for c in ob if c.startswith('program_')]:programs.append({'cluster':c,'program':name[8:],'mean':float(g[name].mean()),'median':float(g[name].median())})
r=pd.concat(ranks);save(r,'unified_descriptive_markers.tsv.gz');save(r.groupby('cluster',sort=False).head(30),'unified_top30_markers.tsv');save(pd.DataFrame(programs),'primary_program_profiles.tsv')
old=ad.read_h5ad(P/'phase9_joint_20260919/objects/joint_CAF_integrated.h5ad',backed='r');shared=old.obs_names.intersection(ob.index);cross=pd.DataFrame({'cell_id':shared,'old_joint_cluster':old.obs.loc[shared,'joint_cluster'].astype(str).to_numpy(),'new_unified_cluster':ob.loc[shared,'unified_cluster'].to_numpy()});save(cross,'old_to_new_retained_cells.tsv');old.file.close()
plt.rcParams.update({'font.size':10})
fig,axs=plt.subplots(2,3,figsize=(17,10));ob['site_display']=ob.site.astype(str).str.lower().replace({'femur':'bone'})
for i,rep in enumerate(['X_umap_unified_uncorrected','X_umap_unified_primary']):
 for j,key in enumerate(['dataset','site_display','unified_cluster']):
  ax=axs[i,j];vals=ob[key].astype(str);cats=sorted(vals.unique(),key=int) if key=='unified_cluster' else sorted(vals.unique())
  for k,v in enumerate(cats):
   m=vals.eq(v).to_numpy();xy=a.obsm[rep][m];ax.scatter(xy[:,0],xy[:,1],s=3,c=[plt.get_cmap('tab20')(k%20)],label='U'+v if key=='unified_cluster' else v,alpha=.8)
  ax.set_title(('Unified input | uncorrected' if i==0 else 'Unified input | Harmony')+' | '+key.replace('_display',''));ax.set_xticks([]);ax.set_yticks([]);ax.legend(fontsize=8,markerscale=3,loc='best')
fig.tight_layout();fig.savefig(O/'figures/01_unified_umap.png',dpi=170);plt.close(fig)
co=changes.groupby('dataset')[['old_joint_cells','unified_strict','osteo_sensitive','competitor_sensitive']].sum();fig,ax=plt.subplots(figsize=(13,6));co.plot.bar(ax=ax);ax.set_ylabel('Candidate cells');ax.set_xlabel('Study');ax.set_title('Previous joint input vs unified gate sets');ax.tick_params(axis='x',rotation=0);fig.tight_layout();fig.savefig(O/'figures/02_gate_changes.png',dpi=180);plt.close(fig)
metrics=pd.read_csv(O/'tables/batch_model_comparison.tsv',sep='\t');mm=metrics[~metrics.analysis.str.startswith('leaveout_')].copy();fig,axs=plt.subplots(1,2,figsize=(15,8));axs[0].barh(mm.analysis,mm.same_batch_neighbors,color='#3a7ca5');axs[0].scatter(mm.composition_expected_same_batch,mm.analysis,color='#c45536',marker='|',s=90,label='Composition expectation');axs[0].invert_yaxis();axs[0].set_xlabel('Same-batch neighbor fraction');axs[0].legend();axs[0].set_title('Batch mixing, interpreted with source composition');axs[1].barh(mm.analysis,mm.ARI_vs_unified_primary,color='#577d5d');axs[1].invert_yaxis();axs[1].set_xlim(0,1.02);axs[1].set_title('Agreement on shared cells');axs[1].set_xlabel('ARI vs unified primary');fig.tight_layout();fig.savefig(O/'figures/03_batch_sensitivity.png',dpi=170);plt.close(fig)
ov=pd.read_csv(O/'tables/model_cluster_overlap.tsv',sep='\t');loo=ov[ov.analysis.str.startswith('leaveout_')];grid=ov[ov.analysis.str.startswith('hvg')];rows=[]
for c in clusters:
 l=loo[(loo.primary_cluster.astype(str)==c)&(loo.primary_cells_shared>=20)];g=grid[grid.primary_cluster.astype(str)==c];rows.append({'cluster':c,'LOO_evaluable_runs':len(l),'LOO_min_Jaccard':l.Jaccard.min(),'LOO_median_Jaccard':l.Jaccard.median(),'HVG_PC_min_Jaccard':g.Jaccard.min(),'HVG_PC_median_Jaccard':g.Jaccard.median()})
st=pd.DataFrame(rows);save(st,'primary_cluster_stability.tsv')
ct=pd.crosstab(ob.unified_cluster,ob.dataset).reindex(clusters);ct=ct.div(ct.sum(1),axis=0);fig,axs=plt.subplots(1,2,figsize=(15,6));ct.plot.bar(stacked=True,ax=axs[0],colormap='tab10');axs[0].set_title('Study composition of unified clusters');axs[0].set_ylabel('Fraction of cells');axs[0].legend(fontsize=8);st.set_index('cluster')[['LOO_min_Jaccard','LOO_median_Jaccard','HVG_PC_min_Jaccard']].plot.bar(ax=axs[1]);axs[1].set_ylim(0,1.03);axs[1].set_title('Cluster stability | shared cells >=20 for leave-out');axs[1].set_ylabel('Best-match Jaccard');fig.tight_layout();fig.savefig(O/'figures/04_cluster_support.png',dpi=170);plt.close(fig)
# Explicit acceptance checks.
assert not ob.independent_doublet_excluded.any();assert ob.ambient_checked.all() and ob.dual_low.all() and ob.reference_QC_pass.all();assert (ob.n_genes>=500).all();assert int(changes.unified_strict.sum())==a.n_obs;assert int(changes.old_joint_cells.sum())==4300;assert all(int((O/'logs'/name).read_text())==0 for name in ['early_ambient.exit_code','unify_gates.exit_code','quality_guard.exit_code','batch_models.exit_code'])
counts_ok=[]
for row in pd.read_csv(O/'resources/object_manifest.tsv',sep='\t').to_dict('records'):
 b=ad.read_h5ad(row['path']);common=a.obs_names.intersection(b.obs_names)
 if len(common):
  bi=b.obs_names.get_indexer(common);ai=a.obs_names.get_indexer(common);vi=b.var_names.get_indexer(a.var_names);assert (vi>=0).all();assert (a.layers['counts'][ai]!=b.layers['counts'][bi][:,vi]).nnz==0;counts_ok.extend(common)
assert set(counts_ok)==set(a.obs_names)
oldman=json.loads((P/'phase9_joint_20260919/MANIFEST_SHA256.json').read_text());assert hashlib.sha256((P/'phase9_joint_20260919/objects/joint_CAF_integrated.h5ad').read_bytes()).hexdigest()==oldman['objects/joint_CAF_integrated.h5ad']
qa={'all_checks_passed':True,'main_cells':a.n_obs,'main_genes':a.n_vars,'main_clusters':len(clusters),'source_matched_raw_counts_cells':len(counts_ok),'old_joint_object_unchanged':True,'primary_dual_CNV_low_and_reference_checked':True,'primary_ambient_checked':True,'reviewed_doublets_in_primary':0,'program_thresholds_remain_operational_not_truth':True}
(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2));(O/'REVIEW_COMPLETE').write_text('UNIFIED_BATCH_GATE_REVIEW_COMPLETE\n');print(json.dumps(qa),flush=True);print(changes.groupby('dataset')[['old_joint_cells','unified_strict','osteo_sensitive','competitor_sensitive','old_removed','new_vs_old']].sum().to_string());print(metrics[['analysis','cells','clusters','same_batch_neighbors','ARI_vs_unified_primary']].to_string(index=False))