from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005';A=B/'phase63_matrix_resistance_20261005'
assert json.loads((O/'RUN_STATUS.json').read_text())['success']
cfg=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text())
for gs in cfg['structure'].values():
 for cs in cfg['contexts'].values():assert set(gs).isdisjoint(cs)
s=pd.read_csv(O/'tables/spatial_gene_and_program_scores.tsv.gz',sep='\t',low_memory=False)
assert len(s)==53309 and s.spot_id.is_unique and s['sample'].nunique()==26 and s.patient_id.nunique()==19
invalid=s[s.library_UMI<=0]
base=pd.read_csv(B/'phase46_spatial_functional_20261001/tables/all_spot_RNA_and_components.tsv.gz',sep='\t',low_memory=False)
assert set(invalid.spot_id)==set(base.loc[base.library_UMI<=0,'spot_id'])
defs={'full':cfg['structure'],'drop_dominant_pair':{k:[g for g in gs if g not in cfg['drop_dominant'][k]] for k,gs in cfg['structure'].items()},'context':cfg['contexts']}
for version,programs in defs.items():
 for name,genes in programs.items():
  cols=['gene__'+g for g in genes if 'gene__'+g in s]
  expected=s[cols].mean(axis=1)
  assert np.allclose(expected,s[version+'__'+name],equal_nan=True)
  assert invalid[version+'__'+name].isna().all()
# Independent rank residual implementation for every raw/depth/mutual sample coefficient.
sample=pd.read_csv(O/'tables/spatial_sample_correlations.tsv',sep='\t',low_memory=False)
def rc(x,y,c):
 arr=np.column_stack([x,y]+c);arr=arr[np.isfinite(arr).all(1)]
 if len(arr)<3:return np.nan
 rr=np.column_stack([rankdata(arr[:,j]) for j in range(arr.shape[1])])
 a=rr[:,0];b=rr[:,1]
 if np.ptp(a)==0 or np.ptp(b)==0:return np.nan
 if c:
  d=np.column_stack([np.ones(len(rr)),rr[:,2:]]);a=a-d@(np.linalg.pinv(d)@a);b=b-d@(np.linalg.pinv(d)@b)
 if np.linalg.norm(a)<1e-9 or np.linalg.norm(b)<1e-9:return np.nan
 return np.corrcoef(a,b)[0,1]
for sm,t in s.groupby('sample'):
 for r in sample[sample['sample']==sm].itertuples():
  c=[] if r.model=='raw' else [np.log1p(t.library_UMI),t.x,t.y]
  if r.model=='mutual_contexts':c += [t['context__'+x] for x in cfg['contexts'] if x!=r.context]
  v=rc(t[r.version+'__'+r.program],t['context__'+r.context],c)
  assert np.isclose(v,r.effect,atol=1e-9,equal_nan=True),(sm,r.program,r.context,r.model,v,r.effect)
# Patient aggregation: equal events then equal sites, no spot-weighted inference.
p=pd.read_csv(O/'tables/spatial_patient_correlations.tsv',sep='\t',low_memory=False);keys=['version','program','context','model']
for scope,d in [('all_patients',sample),('primary_only',sample[sample.site_group=='primary_tumor'])]:
 e=d.groupby(['patient_id','site_group','event_month']+keys,dropna=False).effect.mean()
 v=e.reset_index().groupby(['patient_id','site_group']+keys,dropna=False).effect.mean().reset_index().groupby(['patient_id']+keys,dropna=False).effect.mean().sort_index()
 q=p[p.cohort_scope==scope].set_index(['patient_id']+keys).effect.sort_index()
 assert v.index.equals(q.index) and np.allclose(v,q,equal_nan=True)
 assert q.reset_index().patient_id.nunique()==(19 if scope=='all_patients' else 9)
for filename in ['spatial_patient_statistics.tsv','source_preference_statistics.tsv','discovery_structure_statistics.tsv']:
 d=pd.read_csv(O/'tables'/filename,sep='\t',low_memory=False)
 for fam,t in d.groupby('family'):
  pp=t.P.fillna(1)
  expect=multipletests(pp,method='fdr_bh')[1];ok=t.P.notna()
  assert np.allclose(expect[ok],t.loc[ok,'FDR'])
 assert d[d.N<2].P.isna().all()
# Reproduce inherited full-program subgroup effects.
old=pd.read_csv(A/'tables/posthoc_collagen_structure_statistics.tsv',sep='\t',low_memory=False)
new=pd.read_csv(O/'tables/discovery_structure_statistics.tsv',sep='\t',low_memory=False);new=new[new.version=='full']
keys=['dataset','subgroup_id','matrix','scope','minimum_cells','program']
joined=old.merge(new,on=keys,suffixes=('_old','_new'),validate='one_to_one')
assert len(joined)==len(old)==len(new)
assert np.allclose(joined.mean_effect_old,joined.mean_effect_new,equal_nan=True,atol=1e-10)
assert (joined.N_old==joined.N_new).all()
cells=pd.read_csv(O/'tables/original_cell_assignments_and_confidence.tsv.gz',sep='\t',low_memory=False).fillna('')
oc=pd.read_csv(A/'tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',usecols=list(cells.columns),low_memory=False).fillna('')
pd.testing.assert_frame_equal(cells,oc[cells.columns],check_dtype=False)
# Independent exact sign flips for one complete and one primary-only main spatial endpoint.
def sign_p(v):
 v=np.asarray(v);n=len(v);hits=0
 for i in range(0,2**n,4096):
  ints=np.arange(i,min(i+4096,2**n),dtype=np.uint64)
  signs=1.-2.*((ints[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)
  hits+=sum(np.abs(signs@v)>=abs(v.sum())-1e-11)
 return hits/(2**n)
st=pd.read_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t',low_memory=False)
for scope in ['all_patients','primary_only']:
 q=p[(p.cohort_scope==scope)&(p.version=='full')&(p.program=='fibrillar_I_III_V')&(p.context=='fibro_context')&(p.model=='depth_xy')].effect.dropna()
 r=st[(st.cohort_scope==scope)&(st.version=='full')&(st.program=='fibrillar_I_III_V')&(st.context=='fibro_context')&(st.model=='depth_xy')].iloc[0]
 assert np.isclose(sign_p(q),r.P,atol=1e-12)
# All leave-one-patient-out summaries independently recomputed.
for pf,rf,kk in [('spatial_patient_correlations.tsv','spatial_patient_leave_one_out.tsv',['cohort_scope','version','program','context','model']),('discovery_patient_effects.tsv','discovery_patient_leave_one_out.tsv',['dataset','subgroup_id','matrix','scope','minimum_cells','version','program'])]:
 pp=pd.read_csv(O/'tables'/pf,sep='\t',low_memory=False);rr=pd.read_csv(O/'tables'/rf,sep='\t',low_memory=False).set_index(kk)
 for k,t in pp.groupby(kk):
  vv=t.effect.dropna().to_numpy();r=rr.loc[k]
  if len(vv)>1:
   vals=np.array([np.delete(vv,i).mean() for i in range(len(vv))])
   assert np.isclose(vals.min(),r.loo_min) and np.isclose(vals.max(),r.loo_max)
qa={'passed':True,'all_scores_recalculated':True,'all_sample_correlations_recalculated':True,'patient_aggregation_and_FDR':True,'exact_signflip_independent_spotcheck':True,'inherited_full_effect_rows_reproduced':len(joined),'original8174_assignments_and_confidence_preserved':True,'leave_one_patient_out_recalculated':True,'zero_library_spots_retained':len(invalid),'slices':26,'patient_records':19,'primary_patient_records':9}
(O/'QA/ACCEPTANCE.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa),flush=True)
