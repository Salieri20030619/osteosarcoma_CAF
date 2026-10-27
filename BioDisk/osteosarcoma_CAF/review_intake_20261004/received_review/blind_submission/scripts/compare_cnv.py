import pathlib,json,pandas as pd,numpy as np
P=pathlib.Path(__file__).parent; O=P/'review_output'; dec=json.loads((O/'reviewer_decisions.json').read_text()); s=pd.read_pickle(O/'cell_structure.pkl'); cm=pd.read_csv(P/'copykat_cell_continuous_metrics.tsv.gz',sep='\t'); cm=cm.merge(s[['cell_id','r10','source_id','patient_id']],on='cell_id',validate='one_to_one');cm['label']=cm.r10.map(lambda c:dec[c]['label']);refs=['T','Treg','NK','B','Endo','Lymph','Mast']
ref=cm[cm.label.isin(refs)].groupby('sample_id').agg(ref_cells=('cell_id','size'),ref_clusters=('r10','nunique'),ref_median_abs=('mean_abs_native_signal','median')).reset_index();ref['adequate_reference']=(ref.ref_cells>=30)&(ref.ref_clusters>=2);ref.to_csv(O/'copykat_expression_reference.tsv',sep='\t',index=False)
cc=pd.read_csv(O/'copykat_sample_cluster_metrics.tsv',sep='\t').merge(ref,on='sample_id',how='left');cc['fold_vs_expression_reference']=cc.mean_abs_signal_median/cc.ref_median_abs;cc.loc[~cc.adequate_reference.fillna(False),'fold_vs_expression_reference']=np.nan
cc.loc[cc.copykat_cells<20,'fold_vs_expression_reference']=np.nan
cc.to_csv(O/'copykat_sample_cluster_comparison.tsv',sep='\t',index=False)
patient=cc.groupby(['source_id','cluster_id','patient_id']).agg(copykat_cells=('copykat_cells','sum'),fold=('fold_vs_expression_reference','median'),native_abs=('mean_abs_signal_median','median')).reset_index();patient.to_csv(O/'copykat_patient_comparison.tsv',sep='\t',index=False)
pr=pd.read_csv(O/'copykat_chromosome_profiles.tsv.gz',sep='\t'); pr['label']=pr.cluster_id.map(lambda c:dec[c]['label']);rr=pr[pr.label.isin(refs)].copy();rr['weighted']=rr.native_chromosome_mean*rr.profile_cells;rr=rr.groupby(['sample_id','chromosome']).agg(total=('weighted','sum'),cells=('profile_cells','sum')).reset_index();rr['ref_chromosome_mean']=rr.total/rr.cells
pr=pr.merge(rr[['sample_id','chromosome','ref_chromosome_mean']],on=['sample_id','chromosome'],how='left');pr['ref_centered_signal']=pr.native_chromosome_mean-pr.ref_chromosome_mean
reg=pd.read_csv(O/'regional_CNV_summary.tsv.gz',sep='\t'); cov=pd.read_csv(O/'regional_CNV_coverage.tsv',sep='\t');
# Direction comparison uses non-overlapping chromosome means; event magnitude not interpreted as a common scale.
z=reg[reg.direction.isin(['gain','loss'])].copy();z.loc[z.direction=='loss','value']*=-1
z=z.groupby(['source_id','patient_id','sample_id','cluster_id','method','chromosome','metric'],as_index=False).value.sum().rename(columns={'value':'regional_signed_value'});z=z.merge(pr[['sample_id','cluster_id','chromosome','ref_centered_signal','profile_cells']],on=['sample_id','cluster_id','chromosome'],how='left')
z['eligible_direction_comparison']=(z.ref_centered_signal.abs()>=.02)&(z.regional_signed_value.abs()>=np.where(z.metric=='cell_mean_arm_fraction',.10,10))&(z.profile_cells>=20)
z['direction_agrees']=np.where(z.eligible_direction_comparison,np.sign(z.ref_centered_signal)==np.sign(z.regional_signed_value),np.nan)
z.to_csv(O/'CNV_direction_comparison.tsv.gz',sep='\t',index=False)
cs=[]
for cl,g in patient.groupby('cluster_id'):
 f=g.fold.dropna();rc=cov[cov.cluster_id==cl]; zz=z[(z.cluster_id==cl)&z.eligible_direction_comparison];std=zz[zz.method.isin(['native_HMM','provider_arm'])];
 cs.append({'cluster_id':cl,'copykat_cells':int(g.copykat_cells.sum()),'copykat_patients':len(g),'reference_comparable_patients':len(f),'median_fold':f.median(),'min_fold':f.min(),'max_fold':f.max(),'patients_ge1p5':int((f>=1.5).sum()),'patients_ge2':int((f>=2).sum()),'native_abs_patient_median':g.native_abs.median(),'regional_samples':rc.sample_id.nunique(),'regional_methods':';'.join(sorted(rc.method.unique())),'direction_comparisons':len(std),'direction_agreement':std.direction_agrees.mean()})
cs=pd.DataFrame(cs);cs.to_csv(O/'cnv_cluster_review.tsv',sep='\t',index=False)
print(cs.merge(pd.DataFrame([{'cluster_id':c,'label':d['label']} for c,d in dec.items()]),on='cluster_id').groupby('label').agg(clusters=('cluster_id','size'),median_fold=('median_fold','median'),min_fold=('median_fold','min'),max_fold=('median_fold','max')).to_string())
print('REFERENCE coverage',ref.adequate_reference.sum(),len(ref))
print('COPYKAT source coverage',cm.groupby('source_id').size().to_dict())
print('CNV methods',cov.method.value_counts().to_dict())
print('LARGEST fold',cs.sort_values('median_fold',ascending=False).head(30).to_string(index=False))
(O/'comparison_completed.json').write_text(json.dumps({'status':'complete','clusters_with_copykat':len(cs)}))
