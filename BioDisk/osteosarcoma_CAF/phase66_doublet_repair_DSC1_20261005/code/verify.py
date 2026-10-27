from pathlib import Path
import os,sys,json,hashlib
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
from current_metadata import load_current_cells,apply_current_doublet_evidence
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
new=load_current_cells();old=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
assert new.cell_id.tolist()==old.cell_id.tolist()
protected=['class_id','subgroup_id','lineage_confidence','expression_confidence','functional_interpretation_confidence','CNV_evidence_current','candidate_tier']
for c in protected:pd.testing.assert_series_equal(new[c],old[c],check_names=False)
assert new.doublet_assessed_any.all() and new.doublet_positive_any.sum()==26 and new.legacy_doublet_assessed_any.isna().sum()==5728
# Rebuild union from preserved per-caller fields independently of audit reconciliation.
m=pd.read_csv(O/'tables/doublet_method_details.tsv.gz',sep='\t',low_memory=False).set_index('cell_id').loc[new.cell_id]
assessed=np.zeros(len(m),bool);positive=np.zeros(len(m),bool)
for name in ['input_scDblFinder_class','native_scDblFinder_class']:
 if name in m:
  v=m[name].astype(str).str.lower();ok=v.isin(['singlet','doublet']).to_numpy();assessed|=ok;positive |=ok&v.eq('doublet').to_numpy()
if 'input_predicted_doublet' in m:
 v=m.input_predicted_doublet.astype(str).str.lower();ok=v.isin(['true','false','1','0'])
 # When this field exists only in some sources, NaN status alone does not erase valid native provider flags.
 if 'input_doublet_assessment_status' in m:ok &= ~m.input_doublet_assessment_status.astype(str).str.lower().str.contains('unassess|not_|not |missing|fail|unavail|skip',regex=True)
 assessed|=ok.to_numpy();positive|=ok.to_numpy()&v.isin(['true','1']).to_numpy()
assert np.array_equal(assessed,new.doublet_assessed_any) and np.array_equal(positive,new.doublet_positive_any)
# Prevent recurrence: shuffled IDs align; unknown/duplicate cells fail loudly; no False imputation.
q=apply_current_doublet_evidence(old.sample(frac=1,random_state=19)).set_index('cell_id').loc[new.cell_id]
assert np.array_equal(q.doublet_positive_any,new.doublet_positive_any)
for bad in [pd.DataFrame({'cell_id':['UNKNOWN_CELL']}),pd.DataFrame({'cell_id':[new.cell_id.iloc[0]]*2})]:
 try:apply_current_doublet_evidence(bad)
 except (ValueError,AssertionError):pass
 else:raise AssertionError('Invalid join unexpectedly accepted')
stats=pd.read_csv(O/'tables/program_statistics.tsv.gz',sep='\t');gs=pd.read_csv(O/'tables/gene_statistics.tsv.gz',sep='\t')
for frame in [stats,gs]:
 for fam,t in frame.groupby('family'):
  expected=multipletests(t.P.fillna(1),method='fdr_bh')[1];valid=t.P.notna()
  assert np.allclose(expected[valid],t.loc[valid,'FDR'])
 assert frame[frame.N<2].P.isna().all()
legacy=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/patient_program_statistics.tsv.gz',sep='\t')
keys=['dataset','subgroup_id','matrix','scope','minimum_cells','program']
base=stats[(stats.arm=='keep_all')&(stats.version=='full')&stats.program.isin(legacy.program.unique())]
j=legacy.merge(base,on=keys,suffixes=('_old','_new'),validate='one_to_one')
assert len(j)==len(legacy)==len(base)
assert np.allclose(j.mean_effect_old,j.mean_effect_new,equal_nan=True) and (j.N_old==j.N_new).all()
assert np.allclose(j.P_old,j.P_new,equal_nan=True) and np.allclose(j.FDR_old,j.FDR_new,equal_nan=True)
roster=pd.read_csv(O/'tables/paired_roster.tsv.gz',sep='\t')
n=new.set_index('cell_id')
for r in roster.itertuples():
 t=n[(n.dataset==r.dataset)&(n.donor_record==r.donor_record)&(n.site==r.site)]
 if r.arm=='exclude_model_positive':t=t[~t.doublet_positive_any]
 if r.scope=='Fibro':t=t[t.lineage_stratum=='Fibro']
 assert int((t.subgroup_id==r.subgroup_id).sum())==r.target_cells and int((t.subgroup_id!=r.subgroup_id).sum())==r.reference_cells
 assert min(r.target_cells,r.reference_cells)>=r.minimum_cells
# Independently recalculate affected patient program values from unchanged count arrays.
z=np.load(B/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz');assert list(z['cells'])==list(new.cell_id)
geneidx={x:i for i,x in enumerate(z['genes'])};definitions=json.loads((O/'resources/phase65_protocol.json').read_text())
effects=pd.read_csv(O/'tables/patient_program_effects.tsv.gz',sep='\t')
for state,source in [('F2_C1','GSE270231'),('F2_C3','SCPCP000017')]:
 rr=roster[(roster.arm=='exclude_model_positive')&(roster.dataset==source)&(roster.subgroup_id==state)&(roster.matrix=='corrected')&(roster.scope=='all')&(roster.minimum_cells==1)]
 for donor,t in rr.groupby('donor_record'):
  for program in ['fibrillar_I_III_V','collagen_IV']:
   jj=[geneidx[g] for g in definitions['structure'][program]];vals=[]
   for r in t.itertuples():
    mask=(new.dataset==source)&(new.donor_record==donor)&(new.site==r.site)&~new.doublet_positive_any
    left=np.flatnonzero(mask&new.subgroup_id.eq(state));right=np.flatnonzero(mask&new.subgroup_id.ne(state))
    x=z['corrected'];lib=z['corrected_library']
    dv=np.log2(x[left][:,jj].sum(0)/lib[left].sum()*1e4+.1)-np.log2(x[right][:,jj].sum(0)/lib[right].sum()*1e4+.1)
    vals.append(dv.mean())
   e=effects[(effects.arm=='exclude_model_positive')&(effects.dataset==source)&(effects.subgroup_id==state)&(effects.matrix=='corrected')&(effects.scope=='all')&(effects.minimum_cells==1)&(effects.version=='full')&(effects.program==program)&(effects.donor_record==donor)].effect.iloc[0]
   assert np.isclose(np.mean(vals),e)
# Marker completeness / missingness and actual positive counts.
d=pd.read_csv(O/'tables/desmosome_cell_counts.tsv.gz',sep='\t',low_memory=False,float_precision='round_trip')
r=d[d.matrix=='raw'];assert len(r)==8174 and r.DSC1.notna().sum()==8170 and r.DSC1.sum()==3 and sum(r.DSC1>=1)==3
assert not r[r.DSC1>=1].doublet_positive_any.any()
d[d.DSC1>0].to_csv(O/'tables/DSC1_trace_counts.tsv',sep='\t',index=False,float_format='%.16g')
checks={'passed':True,'all8174_caller_union_rebuilt':True,'fixed_join_order_unknown_duplicate_checks':True,'protected_labels_confidence_unchanged':True,'phase63_baseline_rows_reproduced':len(j),'full_program_and_gene_FDR_verified':True,'every_patient_site_roster_verified':True,'affected_and_main_source_effects_recalculated':True,'DSC1_assayed':8170,'DSC1_raw_detected':3,'retained_sensitivity_cells':8148,'frozen8174_roster_preserved_before_filter':True}
(O/'QA/ACCEPTANCE.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks),flush=True)
