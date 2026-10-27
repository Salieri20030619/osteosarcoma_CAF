from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,h5py,json,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase53_de_novo_reannotation_20261002';B=S/'blind_review';V=P/'review_intake_20261004/received_review';O=P/'phase54_raw_review_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
roster=read(B/'anonymous_cell_cluster_roster.tsv.gz').set_index('cell_id');key=read(S/'unblinding/cell_identifier_key.tsv.gz').set_index('original_cell_id');details=read(V/'blind_submission/cluster_review_details.tsv');labels=details.set_index('cluster_id').lineage_code.to_dict();targets=set(details.loc[details.lineage_code.isin(['Fibro','Mesench'])|details.lineage_code.str.startswith('Mix'),'cluster_id']);refs=['T','Treg','NK','B','Endo','Lymph','Mast']
cm=read(B/'copykat_cell_continuous_metrics.tsv.gz').set_index('cell_id').join(roster[['source_id','patient_id','r10']],validate='one_to_one');assert len(cm)==568754;cm['label']=cm.r10.map(labels)
ref=cm[cm.label.isin(refs)].groupby('sample_id').agg(ref_cells=('r10','size'),ref_clusters=('r10','nunique'),ref_median_abs=('mean_abs_native_signal','median')).reset_index();ref['adequate_reference']=(ref.ref_cells>=30)&(ref.ref_clusters>=2)
cc=cm.groupby(['source_id','r10','sample_id','patient_id']).agg(copykat_cells=('r10','size'),native_abs_median=('mean_abs_native_signal','median')).reset_index().rename(columns={'r10':'cluster_id'}).merge(ref,on='sample_id',how='left')
cc['fold']=cc.native_abs_median/cc.ref_median_abs;cc.loc[~cc.adequate_reference.fillna(False)|(cc.copykat_cells<20),'fold']=np.nan
prior=read(V/'blind_submission/evidence/copykat_sample_cluster_comparison.tsv');check=cc.merge(prior[['sample_id','cluster_id','fold_vs_expression_reference']],on=['sample_id','cluster_id'],validate='one_to_one');assert len(check)==len(cc);assert np.allclose(check.fold,check.fold_vs_expression_reference,equal_nan=True,atol=1e-10)
write(cc,'CNV_sample_recomputed.tsv');pt=cc.groupby(['source_id','cluster_id','patient_id']).agg(copykat_cells=('copykat_cells','sum'),fold=('fold','median')).reset_index();write(pt,'CNV_patient_recomputed.tsv')
# RNA-based genomic calls are kept separately from continuous-expression reference comparisons.
ann=read(S/'tables/all_cell_reannotation_v1.tsv.gz',usecols=['cell_id','copykat_call','infercnv_call','infercnv_source','CNV_evidence_current','identity_interpretation']);ann=ann.rename(columns={'cell_id':'original_cell_id'}).merge(key.reset_index()[['original_cell_id','cell_id']],on='original_cell_id',validate='one_to_one').merge(roster.reset_index(),on='cell_id',validate='one_to_one')
assert len(ann)==876648
calls=ann.groupby(['source_id','r10','patient_id','CNV_evidence_current']).size().rename('cells').reset_index().rename(columns={'r10':'cluster_id'});write(calls,'CNV_patient_native_call_composition.tsv')
write(ann.groupby(['source_id','r10','copykat_call','infercnv_call','CNV_evidence_current'],dropna=False).size().rename('cells').reset_index().rename(columns={'r10':'cluster_id'}),'CNV_cluster_native_call_composition.tsv')
profiles=[];comparisons=[];featurechecks=[];pids=[]
for fi,f in enumerate(sorted((S/'objects/cnv').glob('*_copykat_profiles.h5'))):
 sample=f.name.replace('_copykat_profiles.h5','');bins=read(S/'tables'/f'{sample}_copykat_bins.tsv')
 with h5py.File(f) as h:
  ids=[x.decode() if isinstance(x,bytes) else str(x) for x in h['cell_id'][:]];bn=[x.decode() if isinstance(x,bytes) else str(x) for x in h['bin_id'][:]];x=h['profiles'][:].astype(float)
 assert bn==bins.bin_id.tolist();k=key.loc[ids];obs=roster.loc[k.cell_id].copy();assert obs.sample_id.nunique()==1
 sid=obs.sample_id.iloc[0];source=obs.source_id.iloc[0];patient=obs.patient_id.iloc[0];obs['label']=obs.r10.map(labels);assert np.isfinite(x).all()
 refmask=obs.label.isin(refs).to_numpy();refn=int(refmask.sum());refcl=obs.loc[refmask,'r10'].nunique();adequate=refn>=30 and refcl>=2
 # Rebuild cell-weighted chromosome means directly from the stored cell-by-bin profiles.
 chrs=sorted(set(bins.chr),key=lambda z:(0,int(str(z).replace('chr',''))) if str(z).replace('chr','').isdigit() else (1,str(z)))
 chrX=np.column_stack([x[:,bins.chr.eq(ch).to_numpy()].mean(1) for ch in chrs]);refmean=chrX[refmask].mean(0) if refn else np.full(len(chrs),np.nan)
 refprof=np.median(x[refmask],axis=0) if refn else np.full(x.shape[1],np.nan)
 medians={};nmap={}
 for cl in obs.r10.unique():
  mask=obs.r10.eq(cl).to_numpy();n=int(mask.sum());nmap[cl]=n;medians[cl]=np.median(x[mask],axis=0)-refprof
  mean=chrX[mask].mean(0)
  for ci,ch in enumerate(chrs):profiles.append({'source_id':source,'sample_id':sid,'patient_id':patient,'cluster_id':cl,'chromosome':ch,'cells':n,'native_chromosome_mean':mean[ci],'ref_chromosome_mean':refmean[ci],'ref_centered_signal':mean[ci]-refmean[ci],'adequate_reference':adequate})
 for cl in set(obs.r10)&targets:
  a=medians[cl]
  for other in obs.r10.unique():
   if other==cl:continue
   lineage=labels[other]
   if lineage not in ['Osteo','Osteochondro','Chondro','Fibro','Mural','Myogenic']:continue
   b=medians[other];finite=np.isfinite(a)&np.isfinite(b)
   eligible=adequate and nmap[cl]>=20 and nmap[other]>=20 and finite.sum()>2
   corr=float(np.corrcoef(a[finite],b[finite])[0,1]) if eligible and np.std(a[finite])>0 and np.std(b[finite])>0 else np.nan
   comparisons.append({'source_id':source,'sample_id':sid,'patient_id':patient,'cluster_id':cl,'comparator_cluster':other,'comparator_lineage':lineage,'target_cells':nmap[cl],'comparator_cells':nmap[other],'ref_cells':refn,'ref_clusters':refcl,'eligible_descriptive_comparison':eligible,'residual_bin_Pearson_r':corr,'target_residual_rms':float(np.sqrt(np.mean(a[finite]**2))) if finite.any() else np.nan,'comparator_residual_rms':float(np.sqrt(np.mean(b[finite]**2))) if finite.any() else np.nan,'same_sign_bin_fraction':float(np.mean(np.sign(a[finite])==np.sign(b[finite]))) if eligible else np.nan})
 featurechecks.append({'sample_id':sid,'original_sample':sample,'cells':len(obs),'bins':len(bins),'cell_mapping_and_bin_order_verified':True})
 print(json.dumps({'CNV_samples':fi+1,'sample':sid}),flush=True)
prof=pd.DataFrame(profiles);write(prof,'CNV_chromosome_profiles_recomputed.tsv.gz');write(pd.DataFrame(comparisons),'CNV_same_sample_expression_group_similarity.tsv.gz');write(pd.DataFrame(featurechecks),'CNV_profile_integrity.tsv')
prior=read(V/'blind_submission/evidence/copykat_chromosome_profiles.tsv.gz');ch=prof.merge(prior[['sample_id','cluster_id','chromosome','native_chromosome_mean']],on=['sample_id','cluster_id','chromosome'],suffixes=('_new','_prior'),validate='one_to_one');assert len(ch)==len(prof)
err=float((ch.native_chromosome_mean_new-ch.native_chromosome_mean_prior).abs().max());assert err<1e-6
# Re-evaluate stored region summaries against the rebuilt chromosome-level CopyKAT means.
reg=read(V/'blind_submission/evidence/regional_CNV_summary.tsv.gz');z=reg[reg.direction.isin(['gain','loss'])].copy();z.loc[z.direction.eq('loss'),'value']*=-1
z=z.groupby(['source_id','patient_id','sample_id','cluster_id','method','chromosome','metric'],as_index=False).value.sum().rename(columns={'value':'regional_signed_value'}).merge(prof[['sample_id','cluster_id','chromosome','ref_centered_signal','cells']],on=['sample_id','cluster_id','chromosome'],how='left',validate='many_to_one')
z['eligible_direction_comparison']=(z.ref_centered_signal.abs()>=.02)&(z.regional_signed_value.abs()>=np.where(z.metric.eq('cell_mean_arm_fraction'),.10,10))&(z.cells>=20)
z['direction_agrees']=np.where(z.eligible_direction_comparison,np.sign(z.ref_centered_signal)==np.sign(z.regional_signed_value),np.nan);write(z,'CNV_region_direction_rechecked.tsv.gz')
summary=[]
for cl in sorted(targets):
 q=details[details.cluster_id.eq(cl)].iloc[0];qpt=pt[pt.cluster_id.eq(cl)];val=qpt.fold.dropna();zz=z[(z.cluster_id==cl)&z.eligible_direction_comparison&z.method.isin(['native_HMM','provider_arm'])]
 cp=calls[calls.cluster_id.eq(cl)].groupby('CNV_evidence_current').cells.sum()
 summary.append({'cluster_id':cl,'lineage_code':q.lineage_code,'cells':int(q.cells),'patients_observed':int(q.patients),'copykat_cells':int(qpt.copykat_cells.sum()),'copykat_coverage':qpt.copykat_cells.sum()/q.cells,'comparable_patients':len(val),'median_fold':val.median(),'min_fold':val.min(),'max_fold':val.max(),'patients_ge2_descriptive':int((val>=2).sum()),'direction_records':len(zz),'direction_agreement':zz.direction_agrees.mean(),'native_method_conflict_cells':int(cp.get('method_conflict',0)),'native_concordant_alteration_cells':int(cp.get('concordant_alteration',0)),'native_concordant_no_alteration_cells':int(cp.get('concordant_no_alteration',0))})
write(pd.DataFrame(summary),'CNV_target_summary.tsv')
(O/'CNV_REVIEW_QA.json').write_text(json.dumps({'samples':len(featurechecks),'copykat_cells':len(cm),'all_cell_calls':len(ann),'same_sample_native_amplitude_folds_match_blind_review':True,'chromosome_means_max_difference_from_blind':err,'regional_values':'Reused frozen region summaries; CopyKAT chromosome signals independently rebuilt from cell-by-bin HDF5.','similarity_interpretation':'Descriptive same-method RNA CNV similarity only; no clone identity, DNA confirmation, P value or cell filter.'},indent=2))
(O/'CNV_REVIEW_COMPLETE').write_text('SUCCESS\n');print('CNV_REVIEW_COMPLETE',flush=True)