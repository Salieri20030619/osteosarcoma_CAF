from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import datetime,json,hashlib
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase92_copykat_validation_20261015';O=P/'phase93_expanded_admission_20261016'
for n in ['tables','objects','resources','figures','code','logs']: (O/n).mkdir(parents=True,exist_ok=True)
def dump(f,d):f.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as r:
  for b in iter(lambda:r.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
baseline_path=P/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz'
baseline=pd.read_csv(baseline_path,sep='\t',encoding='utf-8').set_index('cell_id')
e=pd.read_csv(S/'patient_count_verification/cell_count_evidence.tsv.gz',sep='\t',encoding='utf-8').set_index('cell_id')
old=pd.read_csv(P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz',sep='\t',encoding='utf-8').set_index('cell_id')
selected=e[e.reference_validated_reviewed_candidate].copy()
assert len(baseline)==8174 and len(selected)==806 and selected.index.is_unique and not set(selected.index)&set(baseline.index)
assert not set(selected.patient_id)&set(baseline.patient_id.dropna())
extra=[c for c in old.columns if c not in selected.columns]
selected=selected.join(old[extra],validate='one_to_one')
assert selected.expanded_review_retained.all() and not selected.audit_doublet_positive.any()
assert selected.native_copykat_call.eq('diploid').all() and selected.full_scDblFinder_class.eq('singlet').all()
assert selected.provisional_lineage.eq('Fibroblast').all()
assert not selected.manual_cluster_decision.eq('competing_lineage_background').any()
selected['historical_formal_admission']=selected.formal_admission
selected['admitted_to_expanded_resource']=True
selected['confirmed_pure_CAF']=False
selected['site']='Unknown'
selected['admission_scope']='source_stratified_expanded_exploratory'
selected['strict_layer_required_for_entry']=False
selected['new_state_assignment']='pending_reintegration'
selected['HMM_role']='descriptive_only_any_segment_binary_failed_controls'
selected['ambient_evidence_scope']='prior_retained_input_decontX_sensitivity_not_full_input_decontX'
selected.to_csv(O/'tables/admitted_expanded_cells.tsv.gz',sep='\t',encoding='utf-8',index_label='cell_id')
manual=pd.read_csv(P/'audit_expanded_admission_20261015/tables/manual_cluster_adjudication.tsv',sep='\t',encoding='utf-8')
clusters=selected.groupby(['sample','cluster']).agg(cells=('sample','size'),ambient_supported=('ambient_supported_reference_validated_candidate','sum'),legacy_strict=('legacy_strict_supported_reference_validated_candidate','sum'),CTA_context=('CTA_context','sum'),mixed_lineage=('mixed_lineage_flag','sum')).reset_index()
clusters=clusters.merge(manual[['sample','cluster','original_lineage','top15','manual_cluster_decision','manual_reason']],on=['sample','cluster'],validate='one_to_one')
clusters.to_csv(O/'tables/cluster_review.tsv',sep='\t',encoding='utf-8',index=False)
roster=selected.groupby('patient_id').agg(libraries=('sample',lambda x:';'.join(sorted(set(x)))),cells=('sample','size'),ambient_supported=('ambient_supported_reference_validated_candidate','sum'),legacy_strict=('legacy_strict_supported_reference_validated_candidate','sum')).reset_index()
roster['source']='NCC_Zheng2024';roster['site']='Unknown';roster['admission_scope']='expanded_exploratory';roster['Bone_Lung_comparison_eligible']=False
assert len(roster)==6 and roster.cells.sum()==806
roster.to_csv(O/'tables/admitted_patients.tsv',sep='\t',encoding='utf-8',index=False)
registry=baseline.copy();registry['registry_origin']='historical_8174';registry['admitted_to_expanded_resource']=True
new=selected.copy();new['registry_origin']='NCC_addition_806'
# Preserve all old labels; new cells do not inherit old state labels from cluster IDs.
new['source_review_cluster']=new['cluster'];new['cluster']=np.nan
registry=pd.concat([registry,new],axis=0,sort=False)
assert len(registry)==8980 and registry.index.is_unique
registry.to_csv(O/'tables/expanded_registry_8980.tsv.gz',sep='\t',encoding='utf-8',index_label='cell_id')
print('REGISTRY_WRITTEN',len(registry),flush=True)
objects=[];chrom_rows=[];cnvmetrics=[];markers=[]
for sm,g in selected.groupby('sample'):
 print('REVIEW_SAMPLE',sm,flush=True)
 source=P/'phase87_NCC_preprocessing_20261014/inputs'/sm/'counts.h5ad'
 a=ad.read_h5ad(source)
 ids=pd.Index([sm+'_'+str(x) for x in a.obs_names]);assert ids.is_unique
 loc=ids.get_indexer(g.index);assert (loc>=0).all()
 x=a.X[loc].tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all()
 obs=g[['sample','patient_id','source','site','native_copykat_call','full_scDblFinder_class','ambient_supported_reference_validated_candidate','legacy_strict_supported_reference_validated_candidate','admitted_to_expanded_resource','confirmed_pure_CAF','new_state_assignment']].copy()
 b=ad.AnnData(X=x,obs=obs,var=a.var.copy());b.uns['count_semantics']='Unmodified source raw counts; all source genes retained.'
 dest=O/'objects'/(sm+'_admitted_raw_counts.h5ad');b.write_h5ad(dest,compression='gzip')
 read=ad.read_h5ad(dest);assert read.obs_names.equals(g.index) and (read.X!=x).nnz==0
 objects.append({'sample':sm,'path':str(dest),'cells':len(g),'genes':a.n_vars,'sha256':sha(dest),'source':str(source),'raw_count_equality':True})
 if sm=='NCC__BS3':
  ann=pd.read_csv(P/'phase88_expansion_identity_20261014/tables/NCC__BS3_identity.tsv.gz',sep='\t',encoding='utf-8').set_index('cell_id')
  target=ann.index[ann.cluster.astype(str).eq('4')];pos=ids.get_indexer(target);assert (pos>=0).all()
  m=a.X[pos].tocsr();tot=np.asarray(m.sum(axis=1)).ravel();z=m.multiply((1e4/np.maximum(tot,1))[:,None]).tocsr();z.data=np.log1p(z.data)
  means=np.asarray(z.mean(axis=0)).ravel();det=np.asarray((m>0).mean(axis=0)).ravel()
  names=a.var['symbol'].astype('string').fillna('').astype(str)
  q=pd.DataFrame({'gene_id':a.var_names,'symbol':names.values,'mean_log1p_CP10k':means,'detected_fraction':det})
  q=q.sort_values('mean_log1p_CP10k',ascending=False);q.to_csv(O/'tables/BS3_cluster4_descriptive_expression.tsv',sep='\t',encoding='utf-8',index=False)
  panel=['DCN','LUM','COL1A1','COL1A2','COL3A1','COL14A1','PDGFRA','FAP','SFRP2','ASPN','RGS5','PDGFRB','ACTA2','TAGLN','PTPRC','CD3D','LYZ','PECAM1','VWF','EPCAM','MYL1','TTN','ALPL','RUNX2','BGLAP','IBSP','SPP1','IFITM5']
  q[q.symbol.isin(panel)].to_csv(O/'tables/BS3_cluster4_identity_panel.tsv',sep='\t',encoding='utf-8',index=False)
 del a,b,read,x
 c=pd.read_csv(S/'tables'/(sm+'_copykat_full_evidence.tsv.gz'),sep='\t',encoding='utf-8').set_index('cell_id')
 # Whole-chromosome summaries are diagnostics, not additional admission gates.
 controls=c[c.heldout_clean & c.native_copykat_call.isin(['diploid','aneuploid'])]
 ane=c[c.native_copykat_call.eq('aneuploid') & ~c.heldout_clean].sort_index().head(100)
 chosen=pd.concat([c.loc[g.index],controls,ane]);chosen=chosen[~chosen.index.duplicated()]
 toolmap=dict(zip(chosen.tool_id,chosen.index))
 f=S/'copykat'/sm/(sm+'_copykat_final_results_bin_by_cell.txt')
 frame=pd.read_csv(f,sep='\t',encoding='utf-8',usecols=lambda x:x=='chrom' or x in toolmap)
 means=frame.groupby('chrom').mean(numeric_only=True).T.rename(index=toolmap)
 assert set(g.index)<=set(means.index)
 for label,group in [('candidate',g.index),('heldout_normal',controls.index),('aneuploid_context',ane.index)]:
  idx=means.index.intersection(group)
  if len(idx):
   for chrom,val in means.loc[idx].median().items():chrom_rows.append({'sample':sm,'group':label,'chrom':str(chrom),'median_signal':float(val),'cells':len(idx)})
   amp=means.loc[idx].abs().max(axis=1)
   for cell,value in amp.items():cnvmetrics.append({'sample':sm,'cell_id':cell,'group':label,'max_abs_chromosome_mean':float(value)})
pd.DataFrame(chrom_rows).to_csv(O/'tables/chromosome_profile_summary.tsv',sep='\t',encoding='utf-8',index=False)
pd.DataFrame(cnvmetrics).to_csv(O/'tables/chromosome_profile_cell_diagnostics.tsv.gz',sep='\t',encoding='utf-8',index=False)
dump(O/'objects/COUNT_COLLECTION.json',{'baseline':str(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad'),'baseline_cells':8174,'added_objects':objects,'registry':str(O/'tables/expanded_registry_8980.tsv.gz'),'total_cells':8980,'gene_policy':'Separate source matrices retain assayed feature sets. No unassayed gene is imputed as measured zero. Use source-aware gene intersection/availability for downstream analysis.','no_state_transfer':True})
checks=[]
saved=json.loads((P/'phase90_patient_validation_20261015/resources/BASELINE_HASH_CHECK.json').read_text(encoding='utf-8'))['checks']
for path,row in saved.items():
 current=sha(path);checks.append({'path':path,'expected':row['current_sha256'],'actual':current,'pass':current==row['current_sha256']})
for row in json.loads((P/'phase90_patient_validation_20261015/resources/FULL_INPUT_ACCEPTANCE.json').read_text(encoding='utf-8')):
 current=sha(row['source']);checks.append({'path':row['source'],'expected':row['source_sha256'],'actual':current,'pass':current==row['source_sha256']})
assert all(x['pass'] for x in checks)
dump(O/'resources/FINAL_INPUT_INTEGRITY.json',checks)
fig,ax=plt.subplots(figsize=(9,4.8));names=roster.patient_id.str.replace('NCC__','',regex=False)
ax.bar(names,roster.cells,color='#aac9e2',label='Admitted expanded cells');ax.bar(names,roster.ambient_supported,color='#285b80',label='Ambient-supported subset')
for i,n in enumerate(roster.cells):ax.text(i,n+5,str(n),ha='center')
ax.set_ylabel('Cells (not biological replicates)');ax.set_title('NCC expanded addition: 806 cells / 6 source-resolved patients');ax.legend();fig.tight_layout();fig.savefig(O/'figures/expanded_patient_counts.png',dpi=180);plt.close(fig)
q=pd.DataFrame(chrom_rows);q['row']=q['sample']+' | '+q['group'];mat=q.pivot(index='row',columns='chrom',values='median_signal');mat=mat.reindex(columns=sorted(mat.columns,key=lambda x:float(x)))
fig,ax=plt.subplots(figsize=(12,8));im=ax.imshow(mat.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-.15,vmax=.15);ax.set_yticks(range(len(mat)),mat.index,fontsize=8);ax.set_xticks(range(len(mat.columns)),mat.columns);ax.set_xlabel('Chromosome');ax.set_title('Native CopyKAT chromosome summaries (descriptive; no added gate)');fig.colorbar(im,ax=ax,label='Median chromosome mean signal');fig.tight_layout();fig.savefig(O/'figures/chromosome_profiles.png',dpi=160);plt.close(fig)
dump(O/'BUILD_COMPLETE.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_status':'build_complete_final_review_pending','baseline_cells':8174,'added_cells':806,'registry_cells':8980,'NCC_source_patients':6,'ambient_supported_cells':int(roster.ambient_supported.sum()),'legacy_strict_cells':int(roster.legacy_strict.sum()),'integrity_checks_passed':len(checks),'final_review_pending':True})
print('BUILD_COMPLETE',flush=True)
