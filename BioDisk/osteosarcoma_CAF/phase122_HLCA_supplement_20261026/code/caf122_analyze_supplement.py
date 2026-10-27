from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,pandas as pd,numpy as np,anndata as ad
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026'
assert (O/'REFERENCE_HTTP_COMPLETE.json').exists()
ref=ad.read_h5ad(O/'objects/reference_http_program_counts.h5ad');q=ad.read_h5ad(O/'objects/query_program_counts.h5ad');programs=json.loads((O/'resources/SUPPLEMENT_PROGRAMS.json').read_text());common=sorted(set(q.var_names)&set(ref.var_names));q=q[:,common].copy();ref=ref[:,common].copy()
def normalized(a,den):
 x=a.X.astype(float).multiply((10000/np.asarray(den))[:,None]).tocsr();x.data=np.log1p(x.data);return x
qr=normalized(q,q.obs.library_counts);rr=normalized(ref,ref.obs.raw_sum);qs=q.obs.copy();rs=ref.obs.copy();coverage={}
for name,gs in programs.items():
 got=[g for g in gs if g in common];coverage[name]=dict(defined=len(gs),jointly_measured=len(got),genes=got,missing=sorted(set(gs)-set(got)))
 if len(got)/len(gs)<.8:continue
 ix=pd.Index(common).get_indexer(got);qs[name]=np.asarray(qr[:,ix].mean(1)).ravel();rs[name]=np.asarray(rr[:,ix].mean(1)).ravel()
names=[n for n in programs if n in qs];qs['patient_key']=qs.source.astype(str)+'::'+qs.donor.astype(str);rs['patient_key']=rs.study.astype(str)+'::'+rs.donor_id.astype(str)
rg=rs.groupby(['study','patient_key','lung_condition','ann_finest_level'],observed=True);rd=rg[names].mean().reset_index();rd['cells']=rg.size().to_numpy();rd.to_csv(O/'tables/reference_donor_program_scores.tsv',sep='\t',index=False)
nodes=['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2'];qds=[]
for node in nodes:
 g=qs[qs['leaf' if '_s' in node else 'parent'].eq(node)];gb=g.groupby(['source','patient_key'],observed=True);d=gb[names].mean().reset_index();d['cells']=gb.size().to_numpy();d['node_id']=node;qds.append(d)
qd=pd.concat(qds,ignore_index=True);qd.to_csv(O/'tables/query_candidate_patient_program_scores.tsv',sep='\t',index=False)
background=[]
for qmin in [5,10]:
 for rmin in [20,50]:
  for node in nodes:
   a=qd[(qd.node_id==node)&(qd.cells>=qmin)]
   for (condition,typ),b in rd[rd.cells>=rmin].groupby(['lung_condition','ann_finest_level']):
    for name in names:
     if not len(a) or not len(b):continue
     # Average empirical percentile within each study, then equally across studies.
     percentile=[np.mean([np.mean(g[name].to_numpy()<=v) for _,g in b.groupby('study',observed=True)]) for v in a[name]]
     ctr=b.groupby('study',observed=True)[name].mean().mean()
     background.append(dict(node_id=node,query_min_cells=qmin,reference_min_cells=rmin,condition=condition,reference_type=typ,program=name,query_patients=len(a),reference_patients=b.patient_key.nunique(),reference_studies=b.study.nunique(),query_patient_median=float(a[name].median()),reference_equal_study_mean=float(ctr),difference=float(a[name].median()-ctr),median_reference_percentile=float(np.median(percentile)),query_patients_above_ref_mean=int((a[name]>ctr).sum())))
pd.DataFrame(background).to_csv(O/'tables/normal_IPF_background_comparison.tsv',sep='\t',index=False)
tests=[]
for cutoff in [20,50]:
 for (study,typ),g in rd[rd.cells>=cutoff].groupby(['study','ann_finest_level']):
  a=g[g.lung_condition=='IPF'];b=g[g.lung_condition=='Healthy'];shared=set(a.patient_key)&set(b.patient_key)
  a=a[~a.patient_key.isin(shared)];b=b[~b.patient_key.isin(shared)]
  for name in names:
   effect=float(a[name].mean()-b[name].mean()) if len(a) and len(b) else np.nan
   p=float(mannwhitneyu(a[name],b[name],alternative='two-sided',method='auto').pvalue) if min(len(a),len(b))>=3 else np.nan
   tests.append(dict(reference_min_cells=cutoff,study=study,reference_type=typ,program=name,ipf_patients=len(a),healthy_patients=len(b),shared_donors_excluded=len(shared),ipf_minus_healthy=effect,p=p))
tests=pd.DataFrame(tests);tests['fdr']=np.nan
for (cut,study),g in tests.groupby(['reference_min_cells','study']):
 idx=g.index[g.p.notna()]
 if len(idx):tests.loc[idx,'fdr']=multipletests(tests.loc[idx,'p'],method='fdr_bh')[1]
tests['family']='Within reference study and min-cell sensitivity; all assessed programs and stromal types';tests.to_csv(O/'tables/IPF_within_study_patient_tests.tsv',sep='\t',index=False)
ids=pd.read_csv(O/'tables/candidate_identity_patient_summary.tsv',sep='\t');summary=[]
for cutoff in [5,10]:
 for (node,typ),g in ids[ids.cells>=cutoff].groupby(['node_id','identity']):
  summary.append(dict(node_id=node,minimum_candidate_cells=cutoff,identity=typ,patient_records=len(g),equal_patient_fraction=float(g.predicted_fraction.mean()),equal_patient_vote=float(g.mean_vote.mean()),sources=';'.join(sorted(g.source.unique()))))
pd.DataFrame(summary).to_csv(O/'tables/identity_patient_threshold_summary.tsv',sep='\t',index=False)
qs.to_csv(O/'tables/query_cell_program_scores.tsv.gz',sep='\t',index=False)
rs.to_csv(O/'tables/reference_joint_cell_program_scores.tsv.gz',sep='\t',index=False)
(O/'resources/JOINT_PROGRAM_COVERAGE.json').write_text(json.dumps(coverage,indent=2),encoding='utf-8')
rec=dict(status='PASS_COMPUTATION',query_cells=q.n_obs,reference_cells=ref.n_obs,query_candidate_patient_counts={n:{str(k):int(((qd.node_id==n)&(qd.cells>=k)).sum()) for k in [1,5,10]} for n in nodes},reference_patients_by_condition=rs.groupby('lung_condition').patient_key.nunique().to_dict(),reference_studies_by_condition=rs.groupby('lung_condition').study.nunique().to_dict(),joint_genes=len(common),evaluated_programs=names,valid_within_study_tests=int(tests.p.notna().sum()),between_dataset_comparison='Descriptive; no cross-dataset P values',donor_unit='source::donor query; study::donor reference',reference_sampling_cap=100,IPF_signature_source='Sikkema2023 Fig6h CCL2 COL1A1 CTHRC1 MMP19 SERPINE1 HIF1A')
(O/'ANALYSIS_COMPLETE.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec))
